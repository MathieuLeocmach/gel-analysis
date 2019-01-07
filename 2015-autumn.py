from matplotlib.pyplot import *
import sys, os.path
import numpy as np
from colloids import experiment as xp
from colloids.progressbar import ProgressBar
get_ipython().magic(u'matplotlib inline')
get_ipython().system(u'find /media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data -name *_broken_t001.length')
get_ipython().system(u'ls -F --color /media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/1_DiluteGel/172A/1_percolation/*.traj')
get_ipython().system(u'ls -F --color /media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/1_DiluteGel/172A/1_percolation/172A_1206_percolation.*')
pa = '/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data'
ns = [
    '3_DenseGel/163A/1_percolation/1340_percolation',
    '2_MidGel/162B/1_percolation/1745',
    '1_DiluteGel/172A/1_percolation/172A_1206_percolation',
    #'1_DiluteGel/155C/1_percolation/155C_percolation_1645',
    '2_MidGel/170A/1_percolation/170A_1451_percolation',
    '2_MidGel/153A/1_percolation/153A_percolation_1540',
    '3_DenseGel/150A/1_percolation/150A_percolation_1227',
    '3_DenseGel/150A/1_percolation/150A_percolation_1239',
    ]
xs = [
    xp.Experiment(os.path.join(pa, n+'.traj'))
    for n in ns
    ]
hists =[]
for x in xs:
    lengths = [
        np.loadtxt(name, int) 
        for t, name in x.enum('_broken', 'length') 
        if t>0]
    maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
    hists.append(np.array([
        np.histogram(l, np.arange(-1, maxlength+1))[0]
        for l in lengths]))
for his, t0, x in zip(hists, [68,151,18,55,200,104], xs):
    plot(his[t0:,4:].sum(-1)*1./x.get_nb()[t0:-1])
plot(np.arange(len(hists[-1]))+8, hists[-1][:,4:].sum(-1)*1./xs[-1].get_nb()[:-1], gca().lines[-1].get_color())
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
#yscale('log')
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle.pdf')
for his, t0, x in zip(hists, [68,151,18,55,200,104], xs):
    qmax = np.loadtxt(os.path.join(x.path, os.path.splitext(x.trajfile)[0]+'.qmax'))[:,1]
    plot(his[t0:,4:].sum(-1)*1./x.get_nb()[t0:-1]/qmax[t0:len(his)])
qmax = np.loadtxt(os.path.join(xs[-1].path, os.path.splitext(xs[-1].trajfile)[0]+'.qmax'))[:,1]
plot(np.arange(len(hists[-1]))+8, hists[-1][:,4:].sum(-1)*1./xs[-1].get_nb()[:-1]/qmax[:len(hists[-1])], gca().lines[-1].get_color())
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
yscale('log')
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
def smh(h, thr=100):
    """Smooth time dependant histograms by accumulation of at least thr counts. Returns times and counts."""
    dt=0
    ac = 0
    ts = []
    h2 = []
    for i, c in enumerate(h):
        if c > thr:
            if dt>0:
                h2.append(ac*1./dt)
                ts.append(i-dt) #note the starting time of the interval
                ac = 0
                dt = 0
            h2.append(c)
            ts.append(i)
        else:
            ac += c
            dt += 1
            if ac > thr or i+1==len(h):
                h2.append(ac*1./dt)
                ts.append(i+1-dt) #note the starting time of the interval
                ac = 0
                dt = 0
    return np.array(ts), np.array(h2)
qmaxs = []
names = [
    '168A_2004_percolation',
    '172A_1206_percolation',
    '162B_percolation',
    '163A_1340_percolation',
    '162B_1815_ageing'
    ]
imax = 44
for name in names:
    Ss = np.load(os.path.join(
        '/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors',
        name+'_structure_factor.npy'
        ))
    q = np.arange(Ss.shape[1])*Ss[0,0]*8
    qmax.append((Ss[:,4:imax]*q[4:imax]).sum(1)/Ss[:,4:imax].sum(1))
qmaxs = []
names = [
    '168A_2004_percolation',
    '172A_1206_percolation',
    '162B_percolation',
    '163A_1340_percolation',
    '162B_1815_ageing'
    ]
imax = 44
for name in names:
    Ss = np.load(os.path.join(
        '/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors',
        name+'_structure_factor.npy'
        ))
    q = np.arange(Ss.shape[1])*Ss[0,0]*8
    qmaxs.append((Ss[:,4:imax]*q[4:imax]).sum(1)/Ss[:,4:imax].sum(1))
figure(figsize=(6,12))
a = subplot(2,1,1)
for his, t0, x in zip(hists, [68,151,18,55,200,104], xs):
    ts, h = smh(his[t0:,4:].sum(-1))
    plot(ts, h*1./x.get_nb()[t0:-1][ts])
ts, h = smh(hists[-1][:,4:].sum(-1))
plot(ts+8, h*1./xs[-1].get_nb()[:-1][ts], gca().lines[-1].get_color())
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax, t0, x in zip(qmaxs, [68,151,18,55,200,104], xs):
    plot(qmax[t0:])
    
yscale('log');
get_ipython().system(u'ls -F --color /media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors/163A*_structure_factor.npy')
get_ipython().system(u'ls -F --color /media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors/162B*_structure_factor.npy')
get_ipython().system(u'ls -F --color /media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors/172A*_structure_factor.npy')
get_ipython().system(u'ls -F --color /media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors/170A*_structure_factor.npy')
get_ipython().system(u'ls -F --color /media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors/153A*_structure_factor.npy')
get_ipython().system(u'ls -F --color /media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors/150A*_structure_factor.npy')
qmaxs = []
names = [
    '163A_1340_percolation',
    '162B_percolation',
    '172A_1206_percolation',
    '170A_1451_percolation',
    '153A_percolation_1540',
    '150A_percolation_1227',
    '150A_percolation_1239'
    ]
imax = 44
for name in names:
    Ss = np.load(os.path.join(
        '/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors',
        name+'_structure_factor.npy'
        ))
    q = np.arange(Ss.shape[1])*Ss[0,0]*8
    qmaxs.append((Ss[:,4:imax]*q[4:imax]).sum(1)/Ss[:,4:imax].sum(1))
figure(figsize=(6,12))
a = subplot(2,1,1)
for his, t0, x in zip(hists, [68,151,18,55,200,104], xs):
    ts, h = smh(his[t0:,4:].sum(-1))
    plot(ts, h*1./x.get_nb()[t0:-1][ts])
ts, h = smh(hists[-1][:,4:].sum(-1))
plot(ts+8, h*1./xs[-1].get_nb()[:-1][ts], gca().lines[-1].get_color())
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax, t0, x in zip(qmaxs, [68,151,18,55,200,104], xs):
    plot(qmax[t0:])
    
yscale('log');
figure(figsize=(6,12))
t0s = [68,151,18,55,200,104]
a = subplot(2,1,1)
for his, t0, x in zip(hists, t0s, xs):
    ts, h = smh(his[t0:,4:].sum(-1))
    plot(ts, h*1./x.get_nb()[t0:-1][ts])
ts, h = smh(hists[-1][:,4:].sum(-1))
plot(ts+8, h*1./xs[-1].get_nb()[:-1][ts], gca().lines[-1].get_color())
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax, t0, x in zip(qmaxs, t0s, xs):
    plot(qmax[t0:])
plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
yscale('log');
figure(figsize=(6,12))
t0s = [68,151,18,55,200,104]
a = subplot(2,1,1)
for his, t0, x in zip(hists, [68,151,18,55,200,104], xs):
    qmax = np.loadtxt(os.path.join(x.path, os.path.splitext(x.trajfile)[0]+'.qmax'))[:,1]
    plot(his[t0:,4:].sum(-1)*1./x.get_nb()[t0:-1]/qmax[t0:len(his)])
qmax = np.loadtxt(os.path.join(xs[-1].path, os.path.splitext(xs[-1].trajfile)[0]+'.qmax'))[:,1]
plot(np.arange(len(hists[-1]))+8, hists[-1][:,4:].sum(-1)*1./xs[-1].get_nb()[:-1]/qmax[:len(hists[-1])], gca().lines[-1].get_color())
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
plot(0.3*np.arange(1000)**(-4/3.), 'k--')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax, t0, x in zip(qmaxs, t0s, xs):
    plot(qmax[t0:])
plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
yscale('log');
figure(figsize=(6,12))
t0s = [68,151,18,55,200,104]
a = subplot(2,1,1)
for his, t0, x in zip(hists, [68,151,18,55,200,104], xs):
    plot(his[t0:,4:].sum(-1)*1./x.get_nb()[t0:-1])
plot(np.arange(len(hists[-1]))+8, hists[-1][:,4:].sum(-1)*1./xs[-1].get_nb()[:-1], gca().lines[-1].get_color())
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
plot(0.3*np.arange(1000)**(-4/3.), 'k--')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax, t0, x in zip(qmaxs, t0s, xs):
    plot(qmax[t0:])
plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
yscale('log');
figure(figsize=(6,12))
t0s = [75,159,33,55,200,104]
a = subplot(2,1,1)
for his, t0, x in zip(hists, [68,151,18,55,200,104], xs):
    plot(his[t0:,4:].sum(-1)*1./x.get_nb()[t0:-1])
plot(np.arange(len(hists[-1]))+8, hists[-1][:,4:].sum(-1)*1./xs[-1].get_nb()[:-1], gca().lines[-1].get_color())
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
plot(0.3*np.arange(1000)**(-4/3.), 'k--')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax, t0, x in zip(qmaxs, t0s, xs):
    plot(qmax[t0:])
plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
yscale('log');
figure(figsize=(6,12))
t0s = [75,159,33,55,200,104]#[68,151,18,55,200,104]
a = subplot(2,1,1)
for his, t0, x in zip(hists, t0s, xs):
    ts, h = smh(his[t0:,4:].sum(-1))
    plot(ts, h*1./x.get_nb()[t0:-1][ts])
ts, h = smh(hists[-1][:,4:].sum(-1))
plot(ts+8, h*1./xs[-1].get_nb()[:-1][ts], gca().lines[-1].get_color())
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax, t0, x in zip(qmaxs, t0s, xs):
    plot(qmax[t0:])
plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
yscale('log');
figure(figsize=(6,12))
t0s = [68,159,33,55,200,104]#[68,151,18,55,200,104]
a = subplot(2,1,1)
for his, t0, x in zip(hists, t0s, xs):
    ts, h = smh(his[t0:,4:].sum(-1))
    plot(ts, h*1./x.get_nb()[t0:-1][ts])
ts, h = smh(hists[-1][:,4:].sum(-1))
plot(ts+8, h*1./xs[-1].get_nb()[:-1][ts], gca().lines[-1].get_color())
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax, t0, x in zip(qmaxs, t0s, xs):
    plot(qmax[t0:])
plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
yscale('log');
figure(figsize=(6,12))
t0s = [68,151,33,55,200,104]#[68,151,18,55,200,104]
a = subplot(2,1,1)
for his, t0, x in zip(hists, t0s, xs):
    ts, h = smh(his[t0:,4:].sum(-1))
    plot(ts, h*1./x.get_nb()[t0:-1][ts])
ts, h = smh(hists[-1][:,4:].sum(-1))
plot(ts+8, h*1./xs[-1].get_nb()[:-1][ts], gca().lines[-1].get_color())
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax, t0, x in zip(qmaxs, t0s, xs):
    plot(qmax[t0:])
plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
yscale('log');
figure(figsize=(6,12))
t0s = [68,151,33,55,200,104]#[68,151,18,55,200,104]
a = subplot(2,1,1)
for his, t0, x in zip(hists, t0s, xs):
    ts, h = smh(his[t0:,4:].sum(-1))
    plot(ts, h*1./x.get_nb()[t0:-1][ts])
ts, h = smh(hists[-1][:,4:].sum(-1))
plot(ts+8, h*1./xs[-1].get_nb()[:-1][ts], gca().lines[-1].get_color())
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax, t0, x in zip(qmaxs, t0s, xs):
    plot(qmax[t0:])
plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
t0s = [68,151,33,60,200,104]#[68,151,18,55,200,104]
a = subplot(2,1,1)
for his, t0, x in zip(hists, t0s, xs):
    ts, h = smh(his[t0:,4:].sum(-1))
    plot(ts, h*1./x.get_nb()[t0:-1][ts])
ts, h = smh(hists[-1][:,4:].sum(-1))
plot(ts+8, h*1./xs[-1].get_nb()[:-1][ts], gca().lines[-1].get_color())
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax, t0, x in zip(qmaxs, t0s, xs):
    plot(qmax[t0:])
plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
t0s = [68,151,33,65,200,104]#[68,151,18,55,200,104]
a = subplot(2,1,1)
for his, t0, x in zip(hists, t0s, xs):
    ts, h = smh(his[t0:,4:].sum(-1))
    plot(ts, h*1./x.get_nb()[t0:-1][ts])
ts, h = smh(hists[-1][:,4:].sum(-1))
plot(ts+8, h*1./xs[-1].get_nb()[:-1][ts], gca().lines[-1].get_color())
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax, t0, x in zip(qmaxs, t0s, xs):
    plot(qmax[t0:])
plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
t0s = [68,151,33,62,200,104]#[68,151,18,55,200,104]
a = subplot(2,1,1)
for his, t0, x in zip(hists, t0s, xs):
    ts, h = smh(his[t0:,4:].sum(-1))
    plot(ts, h*1./x.get_nb()[t0:-1][ts])
ts, h = smh(hists[-1][:,4:].sum(-1))
plot(ts+8, h*1./xs[-1].get_nb()[:-1][ts], gca().lines[-1].get_color())
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax, t0, x in zip(qmaxs, t0s, xs):
    plot(qmax[t0:])
plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
t0s = [68,151,33,63,200,104]#[68,151,18,55,200,104]
a = subplot(2,1,1)
for his, t0, x in zip(hists, t0s, xs):
    ts, h = smh(his[t0:,4:].sum(-1))
    plot(ts, h*1./x.get_nb()[t0:-1][ts])
ts, h = smh(hists[-1][:,4:].sum(-1))
plot(ts+8, h*1./xs[-1].get_nb()[:-1][ts], gca().lines[-1].get_color())
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax, t0, x in zip(qmaxs, t0s, xs):
    plot(qmax[t0:])
plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
t0s = [68,151,33,64,200,104]#[68,151,18,55,200,104]
a = subplot(2,1,1)
for his, t0, x in zip(hists, t0s, xs):
    ts, h = smh(his[t0:,4:].sum(-1))
    plot(ts, h*1./x.get_nb()[t0:-1][ts])
ts, h = smh(hists[-1][:,4:].sum(-1))
plot(ts+8, h*1./xs[-1].get_nb()[:-1][ts], gca().lines[-1].get_color())
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax, t0, x in zip(qmaxs, t0s, xs):
    plot(qmax[t0:])
plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
yscale('log');
[l.get_color() for l in gca().lines]
qmaxs = []
names = [
    '163A_1340_percolation',
    '162B_percolation',
    '172A_1206_percolation',
    '170A_1451_percolation',
    '153A_percolation_1540',
    '150A_percolation_1227',
    '150A_percolation_1239'
    ]
imax = 44
for name in names:
    Ss = np.load(os.path.join(
        '/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors',
        name+'_structure_factor.npy'
        ))
    q = np.arange(Ss.shape[1])*Ss[0,0]*8
    qmaxs.append((Ss[:,4:imax]*q[4:imax]).sum(1)/Ss[:,4:imax].sum(1))
for qmax in qmaxs[-3:]:
    qmax /=2
figure(figsize=(6,12))
t0s = [68,151,33,64,200,104]#[68,151,18,55,200,104]
a = subplot(2,1,1)
for his, t0, x in zip(hists, t0s, xs):
    ts, h = smh(his[t0:,4:].sum(-1))
    plot(ts, h*1./x.get_nb()[t0:-1][ts])
ts, h = smh(hists[-1][:,4:].sum(-1))
plot(ts+8, h*1./xs[-1].get_nb()[:-1][ts], gca().lines[-1].get_color())
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax, t0, x in zip(qmaxs, t0s, xs):
    plot(qmax[t0:])
plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
yscale('log');
[l.get_color() for l in gca().lines]
2*xs[0].radius, 2*xs[-1].radius
qmaxs = []
names = [
    '163A_1340_percolation',
    '162B_percolation',
    '172A_1206_percolation',
    '170A_1451_percolation',
    '153A_percolation_1540',
    '150A_percolation_1227',
    '150A_percolation_1239'
    ]
imax = 44
for name,x in zip(names, xs):
    Ss = np.load(os.path.join(
        '/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors',
        name+'_structure_factor.npy'
        ))
    q = np.arange(Ss.shape[1])*Ss[0,0]*2*x.radius
    qmaxs.append((Ss[:,4:imax]*q[4:imax]).sum(1)/Ss[:,4:imax].sum(1))
figure(figsize=(6,12))
t0s = [68,151,33,64,200,104]#[68,151,18,55,200,104]
a = subplot(2,1,1)
for his, t0, x in zip(hists, t0s, xs):
    ts, h = smh(his[t0:,4:].sum(-1))
    plot(ts, h*1./x.get_nb()[t0:-1][ts])
ts, h = smh(hists[-1][:,4:].sum(-1))
plot(ts+8, h*1./xs[-1].get_nb()[:-1][ts], gca().lines[-1].get_color())
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax, t0, x in zip(qmaxs, t0s, xs):
    plot(qmax[t0:])
plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
t0s = [68,151,33,64,200,104]#[68,151,18,55,200,104]
a = subplot(2,1,1)
for his, t0, x in zip(hists, t0s, xs):
    ts, h = smh(his[t0:,4:].sum(-1))
    plot(ts, h*1./x.get_nb()[t0:-1][ts])
ts, h = smh(hists[-1][:,4:].sum(-1))
plot(ts+8, h*1./xs[-1].get_nb()[:-1][ts], gca().lines[-1].get_color())
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(0.3*np.arange(1000)**(-1.), 'k--')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax, t0, x in zip(qmaxs, t0s, xs):
    plot(qmax[t0:])
plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
t0s = [68,151,33,64,200,104]#[68,151,18,55,200,104]
a = subplot(2,1,1)
for his, t0, x in zip(hists, t0s, xs):
    ts, h = smh(his[t0:,4:].sum(-1))
    plot(ts, h*1./x.get_nb()[t0:-1][ts])
ts, h = smh(hists[-1][:,4:].sum(-1))
plot(ts+8, h*1./xs[-1].get_nb()[:-1][ts], gca().lines[-1].get_color())
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(0.3/np.log(np.arange(1000)), 'k.')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax, t0, x in zip(qmaxs, t0s, xs):
    plot(qmax[t0:])
plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
t0s = [68,151,33,64,200,104]#[68,151,18,55,200,104]
a = subplot(2,1,1)
for his, t0, x in zip(hists, t0s, xs):
    ts, h = smh(his[t0:,4:].sum(-1))
    plot(ts, h*1./x.get_nb()[t0:-1][ts])
ts, h = smh(hists[-1][:,4:].sum(-1))
plot(ts+8, h*1./xs[-1].get_nb()[:-1][ts], gca().lines[-1].get_color())
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(0.3/np.log(np.arange(1000)), 'k.')
xscale('log')
#ylim(0,0.025)
yscale('log'); #ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax, t0, x in zip(qmaxs, t0s, xs):
    plot(qmax[t0:])
plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
t0s = [68,151,33,64,200,104]#[68,151,18,55,200,104]
a = subplot(2,1,1)
for his, t0, x in zip(hists, t0s, xs):
    ts, h = smh(his[t0:,4:].sum(-1))
    plot(ts, h*1./x.get_nb()[t0:-1][ts])
ts, h = smh(hists[-1][:,4:].sum(-1))
plot(ts+8, h*1./xs[-1].get_nb()[:-1][ts], gca().lines[-1].get_color())
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(1e-3/np.log(np.arange(1000)), 'k.')
xscale('log')
#ylim(0,0.025)
yscale('log'); #ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax, t0, x in zip(qmaxs, t0s, xs):
    plot(qmax[t0:])
plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
t0s = [68,151,33,64,200,104]#[68,151,18,55,200,104]
a = subplot(2,1,1)
for his, t0, x in zip(hists, t0s, xs):
    ts, h = smh(his[t0:,4:].sum(-1))
    plot(ts, h*1./x.get_nb()[t0:-1][ts])
ts, h = smh(hists[-1][:,4:].sum(-1))
plot(ts+8, h*1./xs[-1].get_nb()[:-1][ts], gca().lines[-1].get_color())
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(1e-2/np.log(np.arange(1000)), 'k.')
xscale('log')
#ylim(0,0.025)
yscale('log'); #ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax, t0, x in zip(qmaxs, t0s, xs):
    plot(qmax[t0:])
plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
t0s = [68,151,33,64,200,104]#[68,151,18,55,200,104]
a = subplot(2,1,1)
for his, t0, x in zip(hists, t0s, xs):
    ts, h = smh(his[t0:,4:].sum(-1))
    plot(ts, h*1./x.get_nb()[t0:-1][ts])
ts, h = smh(hists[-1][:,4:].sum(-1))
plot(ts+8, h*1./xs[-1].get_nb()[:-1][ts], gca().lines[-1].get_color())
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(1e-2/np.log(np.arange(1000)), 'k.')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax, t0, x in zip(qmaxs, t0s, xs):
    plot(qmax[t0:])
plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
t0s = [68,151,33,64,200,104]#[68,151,18,55,200,104]
a = subplot(2,1,1)
for his, t0, x in zip(hists, t0s, xs):
    ts, h = smh(his[t0:,4:].sum(-1))
    plot(ts, h*1./x.get_nb()[t0:-1][ts])
ts, h = smh(hists[-1][:,4:].sum(-1))
plot(ts+8, h*1./xs[-1].get_nb()[:-1][ts], gca().lines[-1].get_color())
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(1e-2/np.arange(1000)/np.log(np.arange(1000))**2, 'k.')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax, t0, x in zip(qmaxs, t0s, xs):
    plot(qmax[t0:])
plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
t0s = [68,151,33,64,200,104]#[68,151,18,55,200,104]
a = subplot(2,1,1)
for his, t0, x in zip(hists, t0s, xs):
    ts, h = smh(his[t0:,4:].sum(-1))
    plot(ts, h*1./x.get_nb()[t0:-1][ts])
ts, h = smh(hists[-1][:,4:].sum(-1))
plot(ts+8, h*1./xs[-1].get_nb()[:-1][ts], gca().lines[-1].get_color())
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(np.arange(2,1000), 1e-1/np.arange(2,1000)/np.log(np.arange(2,1000))**2, 'k.')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax, t0, x in zip(qmaxs, t0s, xs):
    plot(qmax[t0:])
plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
t0s = [68,151,33,64,200,104]#[68,151,18,55,200,104]
a = subplot(2,1,1)
for his, t0, x in zip(hists, t0s, xs):
    ts, h = smh(his[t0:,4:].sum(-1))
    plot(ts, h*1./x.get_nb()[t0:-1][ts])
ts, h = smh(hists[-1][:,4:].sum(-1))
plot(ts+8, h*1./xs[-1].get_nb()[:-1][ts], gca().lines[-1].get_color())
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(np.arange(2,1000), 1./np.arange(2,1000)/np.log(np.arange(2,1000))**2, 'k.')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax, t0, x in zip(qmaxs, t0s, xs):
    plot(qmax[t0:])
plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
t0s = [68,151,33,64,200,104]#[68,151,18,55,200,104]
a = subplot(2,1,1)
for his, t0, x in zip(hists, t0s, xs):
    ts, h = smh(his[t0:,4:].sum(-1))
    plot(ts, h*1./x.get_nb()[t0:-1][ts])
ts, h = smh(hists[-1][:,4:].sum(-1))
plot(ts+8, h*1./xs[-1].get_nb()[:-1][ts], gca().lines[-1].get_color())
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(np.arange(2,1000), 1./np.arange(2,1000)/np.log(np.arange(2,1000)/10.)**2, 'k.')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax, t0, x in zip(qmaxs, t0s, xs):
    plot(qmax[t0:])
plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
t0s = [68,151,33,64,200,104]#[68,151,18,55,200,104]
a = subplot(2,1,1)
for his, t0, x in zip(hists, t0s, xs):
    ts, h = smh(his[t0:,4:].sum(-1))
    plot(ts, h*1./x.get_nb()[t0:-1][ts])
ts, h = smh(hists[-1][:,4:].sum(-1))
plot(ts+8, h*1./xs[-1].get_nb()[:-1][ts], gca().lines[-1].get_color())
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(np.arange(2,1000), 1./np.arange(2,1000)/np.log(np.arange(2,1000)/0.10.)**2, 'k.')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax, t0, x in zip(qmaxs, t0s, xs):
    plot(qmax[t0:])
plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
t0s = [68,151,33,64,200,104]#[68,151,18,55,200,104]
a = subplot(2,1,1)
for his, t0, x in zip(hists, t0s, xs):
    ts, h = smh(his[t0:,4:].sum(-1))
    plot(ts, h*1./x.get_nb()[t0:-1][ts])
ts, h = smh(hists[-1][:,4:].sum(-1))
plot(ts+8, h*1./xs[-1].get_nb()[:-1][ts], gca().lines[-1].get_color())
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(np.arange(2,1000), 1./np.arange(2,1000)/np.log(np.arange(2,1000)/0.10)**2, 'k.')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax, t0, x in zip(qmaxs, t0s, xs):
    plot(qmax[t0:])
plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
t0s = [68,151,33,64,200,104]#[68,151,18,55,200,104]
a = subplot(2,1,1)
for his, t0, x in zip(hists, t0s, xs):
    ts, h = smh(his[t0:,4:].sum(-1))
    plot(ts, h*1./x.get_nb()[t0:-1][ts])
ts, h = smh(hists[-1][:,4:].sum(-1))
plot(ts+8, h*1./xs[-1].get_nb()[:-1][ts], gca().lines[-1].get_color())
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(np.arange(2,1000), 1./np.arange(2,1000)/np.log(np.arange(2,1000))**2, 'k.')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax, t0, x in zip(qmaxs, t0s, xs):
    plot(qmax[t0:])
plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
t0s = [68,151,33,64,200,104]#[68,151,18,55,200,104]
a = subplot(2,1,1)
for his, t0, x in zip(hists, t0s, xs):
    ts, h = smh(his[t0:,4:].sum(-1))
    plot(ts, h*1./x.get_nb()[t0:-1][ts])
ts, h = smh(hists[-1][:,4:].sum(-1))
plot(ts+8, h*1./xs[-1].get_nb()[:-1][ts], gca().lines[-1].get_color())
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(1e-2/np.log(np.arange(1000)), 'k.')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax, t0, x in zip(qmaxs, t0s, xs):
    plot(qmax[t0:])
plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
t0s = [68,151,33,64,200,104]#[68,151,18,55,200,104]
a = subplot(2,1,1)
for his, t0, x in zip(hists, t0s, xs):
    ts, h = smh(his[t0:,4:].sum(-1))
    plot(ts, h*1./x.get_nb()[t0:-1][ts])
ts, h = smh(hists[-1][:,4:].sum(-1))
plot(ts+8, h*1./xs[-1].get_nb()[:-1][ts], gca().lines[-1].get_color())
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax, t0, x in zip(qmaxs, t0s, xs):
    plot(qmax[t0:])
plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
yscale('log');
[l.get_color() for l in gca().lines]
pa = '/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data'
ns = [
    '3_DenseGel/163A/1_percolation/1340_percolation',
    '2_MidGel/162B/1_percolation/1745',
    '1_DiluteGel/172A/1_percolation/172A_1206_percolation',
    #'1_DiluteGel/155C/1_percolation/155C_percolation_1645',
    '2_MidGel/170A/1_percolation/170A_1451_percolation',
    '2_MidGel/153A/1_percolation/153A_percolation_1540',
    '3_DenseGel/150A/1_percolation/150A_percolation_1227',
    '3_DenseGel/150A/1_percolation/150A_percolation_1239',
    '2_MidGel/162B/2_ageing/1815_ageing'
    ]
xs = [
    xp.Experiment(os.path.join(pa, n+'.traj'))
    for n in ns
    ]
hists =[]
for x in xs:
    lengths = [
        np.loadtxt(name, int) 
        for t, name in x.enum('_broken', 'length') 
        if t>0]
    maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
    hists.append(np.array([
        np.histogram(l, np.arange(-1, maxlength+1))[0]
        for l in lengths]))
for his, t0, x in zip(hists, [68,151,18,55,200,104], xs):
    plot(his[t0:,4:].sum(-1)*1./x.get_nb()[t0:-1])
    
plot(np.arange(len(hists[-2]))+8, hists[-2][:,4:].sum(-1)*1./xs[-2].get_nb()[:-1], gca().lines[-1].get_color())
plot(np.arange(len(hists))*3+35*6-159, hists[-1][:,4:].sum(-1)*1./xs[-1].get_nb()[:-1], 'g')

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
#yscale('log')
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle.pdf')
for his, t0, x in zip(hists, [68,151,18,55,200,104], xs):
    plot(his[t0:,4:].sum(-1)*1./x.get_nb()[t0:-1])
    
plot(np.arange(len(hists[-2]))+8, hists[-2][:,4:].sum(-1)*1./xs[-2].get_nb()[:-1], gca().lines[-1].get_color())
plot(np.arange(len(hists[-1]))*3+35*6-159, hists[-1][:,4:].sum(-1)*1./xs[-1].get_nb()[:-1], 'g')

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
#yscale('log')
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle.pdf')
t0s = [75,159,33,55,200,104,0,0]
rates = [
    np.column_stack((
        np.arange(len(his)-t0),
        his[t0:,4:].sum(-1)*1./x.get_nb()[t0:-1]
        )) 
    for his, t0, x in zip(hists, t0s, xs)
    ]

#150A
rates[5] = np.vstack((rates[5], rates[6]+[8,0]))
del rates[6]
#162B
rates[1] = np.vstack((rates[1], rates[5]*[3,0]+[35*6-159,0]))
del rates[5]
for r in rates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
#yscale('log')
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle.pdf')
t0s = [68,151,33,64,200,104,0,0]
rates = [
    np.column_stack((
        np.arange(len(his)-t0),
        his[t0:,4:].sum(-1)*1./x.get_nb()[t0:-1]
        )) 
    for his, t0, x in zip(hists, t0s, xs)
    ]

#150A
rates[5] = np.vstack((rates[5], rates[6]+[8,0]))
del rates[6]
#162B
rates[1] = np.vstack((rates[1], rates[5]*[3,0]+[35*6-159,0]))
del rates[5]
for r in rates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
#yscale('log')
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle.pdf')
t0s = [68,151,33,64,200,104,0,0]
rates = [
    np.column_stack((
        np.arange(len(his)-t0),
        his[t0:,4:].sum(-1)*1./x.get_nb()[t0:-1]
        )) 
    for his, t0, x in zip(hists, t0s, xs)
    ]

#150A
rates[5] = np.vstack((rates[5], rates[6]+[8,0]))
del rates[6]
#162B
rates[1] = np.vstack((rates[1], rates[5]*[3,1]+[35*6-159,0]))
del rates[5]
for r in rates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
#yscale('log')
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle.pdf')
t0s = [68,151,33,64,200,104,0,0]
rates = [
    np.column_stack((
        np.arange(len(his)-t0),
        his[t0:,4:].sum(-1)*1./x.get_nb()[t0:-1]
        )) 
    for his, t0, x in zip(hists, t0s, xs)
    ]

#162B
rates[1] = np.vstack((rates[1], rates[-1]*[3,1]+[35*6-159,0]))
del rates[-1]
#150A
rates[5] = np.vstack((rates[5], rates[6]+[8,0]))
del rates[6]
for r in rates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
#yscale('log')
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle.pdf')
for qmax, t0, x in zip(qmaxs, t0s, xs):
    plot(qmax[t0:])
plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
xscale('log');
yscale('log');
xlabel(r'$t/\tau_B$')
ylabel(r'$q_\mathrm{max}\sigma$')
get_ipython().system(u'ls -F --color /media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors/162B*_structure_factor.npy')
qmaxs = []
names = [
    '163A_1340_percolation',
    '162B_percolation',
    '172A_1206_percolation',
    '170A_1451_percolation',
    '153A_percolation_1540',
    '150A_percolation_1227',
    '150A_percolation_1239',
    '162B_ageing',
    ]
imax = 44
for name, t0, x in zip(names, t0s, xs):
    Ss = np.load(os.path.join(
        '/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors',
        name+'_structure_factor.npy'
        ))
    q = np.arange(Ss.shape[1])*Ss[0,0]*2*x.radius
    qmaxs.append(
        np.column_stack((
            np.arange(len(Ss)-t0),
            (Ss[:,4:imax]*q[4:imax]).sum(1)/Ss[:,4:imax].sum(1)
            ))
        )
#162B
qmaxs[1] = np.vstack((qmaxs[1], qmaxs[-1]*[3,1]+[35*6-159,0]))
del qmaxs[-1]
#150A
qmaxs[5] = np.vstack((qmaxs[5], qmaxs[6]+[8,0]))
del qmaxs[6]
qmaxs = []
names = [
    '163A_1340_percolation',
    '162B_percolation',
    '172A_1206_percolation',
    '170A_1451_percolation',
    '153A_percolation_1540',
    '150A_percolation_1227',
    '150A_percolation_1239',
    '162B_1815_ageing',
    ]
imax = 44
for name, t0, x in zip(names, t0s, xs):
    Ss = np.load(os.path.join(
        '/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors',
        name+'_structure_factor.npy'
        ))
    q = np.arange(Ss.shape[1])*Ss[0,0]*2*x.radius
    qmaxs.append(
        np.column_stack((
            np.arange(len(Ss)-t0),
            (Ss[:,4:imax]*q[4:imax]).sum(1)/Ss[:,4:imax].sum(1)
            ))
        )
#162B
qmaxs[1] = np.vstack((qmaxs[1], qmaxs[-1]*[3,1]+[35*6-159,0]))
del qmaxs[-1]
#150A
qmaxs[5] = np.vstack((qmaxs[5], qmaxs[6]+[8,0]))
del qmaxs[6]
name
qmaxs = []
names = [
    '163A_1340_percolation',
    '162B_percolation',
    '172A_1206_percolation',
    '170A_1451_percolation',
    '153A_percolation_1540',
    '150A_percolation_1227',
    '150A_percolation_1239',
    '162B_1815_ageing',
    ]
imax = 44
for name, t0, x in zip(names, t0s, xs):
    Ss = np.load(os.path.join(
        '/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors',
        name+'_structure_factor.npy'
        ))
    q = np.arange(Ss.shape[1])*Ss[0,0]*2*x.radius
    qmaxs.append(
        np.column_stack((
            np.arange(len(Ss)-t0),
            (Ss[t0:,4:imax]*q[4:imax]).sum(1)/Ss[t0:,4:imax].sum(1)
            ))
        )
#162B
qmaxs[1] = np.vstack((qmaxs[1], qmaxs[-1]*[3,1]+[35*6-159,0]))
del qmaxs[-1]
#150A
qmaxs[5] = np.vstack((qmaxs[5], qmaxs[6]+[8,0]))
del qmaxs[6]
for qmax, t0, x in zip(qmaxs, t0s, xs):
    plot(qmax[t0:])
plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
xscale('log');
yscale('log');
xlabel(r'$t/\tau_B$')
ylabel(r'$q_\mathrm{max}\sigma$')
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
xscale('log');
yscale('log');
xlabel(r'$t/\tau_B$')
ylabel(r'$q_\mathrm{max}\sigma$')
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
#plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
xscale('log');
yscale('log');
xlabel(r'$t/\tau_B$')
ylabel(r'$q_\mathrm{max}\sigma$')
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in rates:
    plot(r[:,0], r[:,1])
    
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
plot(0.3*np.arange(1000)**(-4/3.), 'k--')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
def smh(h, thr=100):
    """Smooth time dependant histograms by accumulation of at least thr counts. Returns times and counts."""
    dt=0
    ac = 0
    ts = []
    h2 = []
    for i, c in enumerate(h):
        if c > thr:
            if dt>0:
                h2.append(ac*1./dt)
                ts.append(i-dt) #note the starting time of the interval
                ac = 0
                dt = 0
            h2.append(c)
            ts.append(i)
        else:
            ac += c
            dt += 1
            if ac > thr or i+1==len(h):
                h2.append(ac*1./dt)
                ts.append(i+1-dt) #note the starting time of the interval
                ac = 0
                dt = 0
    return np.array(ts), np.array(h2)

smrates = [
    np.column_stack(smh(his[t0:,4:].sum(-1))) * np.column_stack((np.ones(len(his)-t0), 1./x.get_nb()[t0:-1]))
    for his, t0, x in zip(hists, t0s, xs)
    ]

#162B
smrates[1] = np.vstack((smrates[1], smrates[-1]*[3,1]+[35*6-159,0]))
del smrates[-1]
#150A
smrates[5] = np.vstack((smrates[5], smrates[6]+[8,0]))
del smrates[6]
def smh(h, thr=100):
    """Smooth time dependant histograms by accumulation of at least thr counts. Returns times and counts."""
    dt=0
    ac = 0
    ts = []
    h2 = []
    for i, c in enumerate(h):
        if c > thr:
            if dt>0:
                h2.append(ac*1./dt)
                ts.append(i-dt) #note the starting time of the interval
                ac = 0
                dt = 0
            h2.append(c)
            ts.append(i)
        else:
            ac += c
            dt += 1
            if ac > thr or i+1==len(h):
                h2.append(ac*1./dt)
                ts.append(i+1-dt) #note the starting time of the interval
                ac = 0
                dt = 0
    return np.array(ts), np.array(h2)

smrates = [
    np.column_stack(smh(his[t0:,4:].sum(-1)))
    for his, t0, x in zip(hists, t0s, xs)
    ]
for r in smrates:
    r[:,1] /= x.get_nb()[t0:][r[:,0]]

#162B
smrates[1] = np.vstack((smrates[1], smrates[-1]*[3,1]+[35*6-159,0]))
del smrates[-1]
#150A
smrates[5] = np.vstack((smrates[5], smrates[6]+[8,0]))
del smrates[6]
def smh(h, thr=100):
    """Smooth time dependant histograms by accumulation of at least thr counts. Returns times and counts."""
    dt=0
    ac = 0
    ts = []
    h2 = []
    for i, c in enumerate(h):
        if c > thr:
            if dt>0:
                h2.append(ac*1./dt)
                ts.append(i-dt) #note the starting time of the interval
                ac = 0
                dt = 0
            h2.append(c)
            ts.append(i)
        else:
            ac += c
            dt += 1
            if ac > thr or i+1==len(h):
                h2.append(ac*1./dt)
                ts.append(i+1-dt) #note the starting time of the interval
                ac = 0
                dt = 0
    return np.array(ts), np.array(h2)

smrates = [
    np.column_stack(smh(his[t0:,4:].sum(-1)))
    for his, t0, x in zip(hists, t0s, xs)
    ]
for r in smrates:
    r[:,1] /= x.get_nb()[t0:][r[:,0].astype(int)]

#162B
smrates[1] = np.vstack((smrates[1], smrates[-1]*[3,1]+[35*6-159,0]))
del smrates[-1]
#150A
smrates[5] = np.vstack((smrates[5], smrates[6]+[8,0]))
del smrates[6]
def smh(h, thr=100):
    """Smooth time dependant histograms by accumulation of at least thr counts. Returns times and counts."""
    dt=0
    ac = 0
    ts = []
    h2 = []
    for i, c in enumerate(h):
        if c > thr:
            if dt>0:
                h2.append(ac*1./dt)
                ts.append(i-dt) #note the starting time of the interval
                ac = 0
                dt = 0
            h2.append(c)
            ts.append(i)
        else:
            ac += c
            dt += 1
            if ac > thr or i+1==len(h):
                h2.append(ac*1./dt)
                ts.append(i+1-dt) #note the starting time of the interval
                ac = 0
                dt = 0
    return np.array(ts), np.array(h2)

smrates = [
    np.column_stack(smh(his[t0:,4:].sum(-1)))
    for his, t0, x in zip(hists, t0s, xs)
    ]
for r in smrates:
    r[:,1] /= x.get_nb()[t0-1:][r[:,0].astype(int)]

#162B
smrates[1] = np.vstack((smrates[1], smrates[-1]*[3,1]+[35*6-159,0]))
del smrates[-1]
#150A
smrates[5] = np.vstack((smrates[5], smrates[6]+[8,0]))
del smrates[6]
r
r[:,0].astype(int)
r[:,0].astype(int)<x.size
def smh(h, thr=100):
    """Smooth time dependant histograms by accumulation of at least thr counts. Returns times and counts."""
    dt=0
    ac = 0
    ts = []
    h2 = []
    for i, c in enumerate(h):
        if c > thr:
            if dt>0:
                h2.append(ac*1./dt)
                ts.append(i-dt) #note the starting time of the interval
                ac = 0
                dt = 0
            h2.append(c)
            ts.append(i)
        else:
            ac += c
            dt += 1
            if ac > thr or i+1==len(h):
                h2.append(ac*1./dt)
                ts.append(i+1-dt) #note the starting time of the interval
                ac = 0
                dt = 0
    return np.array(ts), np.array(h2)

smrates = [
    np.column_stack(smh(his[t0:,4:].sum(-1)))
    for his, t0 in zip(hists, t0s)
    ]
for r, x in zip(smrates, xs):
    r[:,1] /= x.get_nb()[t0:][r[:,0].astype(int)]

#162B
smrates[1] = np.vstack((smrates[1], smrates[-1]*[3,1]+[35*6-159,0]))
del smrates[-1]
#150A
smrates[5] = np.vstack((smrates[5], smrates[6]+[8,0]))
del smrates[6]
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
[l.get_color() for l in gca().lines]
def smh(h, thr=100):
    """Smooth time dependant histograms by accumulation of at least thr counts. Returns times and counts."""
    dt=0
    ac = 0
    ts = []
    h2 = []
    for i, c in enumerate(h):
        if c > thr:
            if dt>0:
                h2.append(ac*1./dt)
                ts.append(i-dt) #note the starting time of the interval
                ac = 0
                dt = 0
            h2.append(c)
            ts.append(i)
        else:
            ac += c
            dt += 1
            if ac > thr or i+1==len(h):
                h2.append(ac*1./dt)
                ts.append(i+1-dt) #note the starting time of the interval
                ac = 0
                dt = 0
    return np.array(ts), np.array(h2)

smrates = [
    np.column_stack(smh(his[t0:,4:].sum(-1)))
    for his, t0 in zip(hists, t0s)
    ]
for r, x in zip(smrates, xs):
    r[:,1] /= x.get_nb()[t0:][r[:,0].astype(int)]

#162B
smrates[1] = np.vstack((smrates[1], smrates[-1]*[1,1]+[35*6-159,0]))
del smrates[-1]
#150A
smrates[5] = np.vstack((smrates[5], smrates[6]+[8,0]))
del smrates[6]
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
[l.get_color() for l in gca().lines]
qmaxs = []
names = [
    '163A_1340_percolation',
    '162B_percolation',
    '172A_1206_percolation',
    '170A_1451_percolation',
    '153A_percolation_1540',
    '150A_percolation_1227',
    '150A_percolation_1239',
    '162B_1815_ageing',
    ]
imax = 44
for name, t0, x in zip(names, t0s, xs):
    Ss = np.load(os.path.join(
        '/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors',
        name+'_structure_factor.npy'
        ))
    q = np.arange(Ss.shape[1])*Ss[0,0]*2*x.radius
    qmaxs.append(
        np.column_stack((
            np.arange(len(Ss)-t0),
            (Ss[t0:,4:imax]*q[4:imax]).sum(1)/Ss[t0:,4:imax].sum(1)
            ))
        )
#162B
qmaxs[1] = np.vstack((qmaxs[1], qmaxs[-1]*[1,1]+[35*6-159,0]))
del qmaxs[-1]
#150A
qmaxs[5] = np.vstack((qmaxs[5], qmaxs[6]+[8,0]))
del qmaxs[6]
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
#plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
xscale('log');
yscale('log');
xlabel(r'$t/\tau_B$')
ylabel(r'$q_\mathrm{max}\sigma$')
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in rates:
    plot(r[:,0], r[:,1])
    
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
plot(0.3*np.arange(1000)**(-4/3.), 'k--')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
t0s = [68,151,33,64,200,104,0,0]
rates = [
    np.column_stack((
        np.arange(len(his)-t0),
        his[t0:,4:].sum(-1)*1./x.get_nb()[t0:-1]
        )) 
    for his, t0, x in zip(hists, t0s, xs)
    ]

#162B
rates[1] = np.vstack((rates[1], rates[-1]*[1,1]+[35*6-159,0]))
del rates[-1]
#150A
rates[5] = np.vstack((rates[5], rates[6]+[8,0]))
del rates[6]
for r in rates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
#yscale('log')
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle.pdf')
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in rates:
    plot(r[:,0], r[:,1])
    
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
plot(0.3*np.arange(1000)**(-4/3.), 'k--')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(0.3*np.arange(1000)**(-0.5), 'k-.')
plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(1e-2*np.arange(1000)**(-0.5), 'k-.')
plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(1e-2*np.arange(1000)**(-1/4.), 'k-.')
plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(1e-2*np.arange(1000)**(-1/3.), 'k-.')
plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(1e-2*np.arange(1000)**(-1/4.), 'k-.')
plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(0.5*(qmax[1:,0]+qmax[:-1,0]), np.diff(qmax[:,1])/np.diff(qmax[:,0]))
yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], gaussian_filter1d(qmax[:,1], 1, order=1))
yscale('log');
[l.get_color() for l in gca().lines]
from scipy.ndimage import gaussian_filter1d
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], gaussian_filter1d(qmax[:,1], 1, order=1))
yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], gaussian_filter1d(qmax[:,1], 2, order=1))
yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], gaussian_filter1d(qmax[:,1], 8, order=1))
yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], gaussian_filter1d(qmax[:,1], 8, order=1))
#yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], -gaussian_filter1d(qmax[:,1], 8, order=1))
#yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], -gaussian_filter1d(qmax[:,1], 8, order=1))
yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], -gaussian_filter1d(qmax[:,1], 1, order=1))
yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a, sharey=b)
for qmax in qmaxs:
    plot(qmax[:,0], -gaussian_filter1d(qmax[:,1], 1, order=1))
#yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a, sharey=a)
for qmax in qmaxs:
    plot(qmax[:,0], -gaussian_filter1d(qmax[:,1], 1, order=1))
#yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a, sharey=a)
for qmax in qmaxs:
    plot(qmax[:,0], -gaussian_filter1d(qmax[:,1], 1, order=1))
#yscale('log');
ylim(1e-3, 4e-2)
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], np.cumsum(r[:,1]))

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a, sharey=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
#yscale('log');
ylim(1e-3, 4e-2)
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], np.cumsum(r[:,1]))

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); #ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
ylim(1e-3, 4e-2)
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], 1./np.cumsum(r[:,1]))

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); #ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
ylim(1e-3, 4e-2)
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], 1./np.cumsum(r[:,1]))

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-1/3.), 'k--')
#plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
#plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); #ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
ylim(1e-3, 4e-2)
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], 1./np.cumsum(r[:,1]))

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-1/3.), 'k--')
#plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
#plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); #ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], np.cumsum(r[::-1,1])[::-1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-1/3.), 'k--')
#plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
#plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); #ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], np.cumsum(r[::-1,1])[::-1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-1/3.), 'k--')
#plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
#plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); #ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], np.cumsum(r[:,1])*1./r[:,1].sum())

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-1/3.), 'k--')
#plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
#plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); #ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], np.cumsum(r[:,1])-r[:,1].sum())

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-1/3.), 'k--')
#plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
#plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); #ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], r[:,1].sum()-np.cumsum(r[:,1]))

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-1/3.), 'k--')
#plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
#plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); #ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
a = subplot(3,1,1)
for r in smrates:
    plot(r[:,0], r[:,1])
subplot(3,1,2, sharex=a)
for r in rates:
    plot(r[:,0], r[:,1].sum()-np.cumsum(r[:,1]))

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-1/3.), 'k--')
#plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
#plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); #ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(3,1,3, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
a = subplot(3,1,1)
yscale('log');
for r in smrates:
    plot(r[:,0], r[:,1])
plot(0.3*np.arange(1000)**(-4/3.), 'k--')

subplot(3,1,2, sharex=a)
for r in rates:
    plot(r[:,0], r[:,1].sum()-np.cumsum(r[:,1]))

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-1/3.), 'k--')
#plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
#plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); #ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(3,1,3, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
a = subplot(3,1,1)
xscale('log'); yscale('log');
for r in smrates:
    plot(r[:,0], r[:,1])
plot(0.3*np.arange(1000)**(-4/3.), 'k--')

subplot(3,1,2, sharex=a)
for r in rates:
    plot(r[:,0], r[:,1].sum()-np.cumsum(r[:,1]))

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-1/3.), 'k--')
#plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
#plot(1e-2/np.log(np.arange(1000)), 'k:')

#ylim(0,0.025)
yscale('log'); #ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(3,1,3, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
a = subplot(3,1,1)
xscale('log'); yscale('log'); ylim(1e-3,0.025)
for r in smrates:
    plot(r[:,0], r[:,1])
plot(0.3*np.arange(1000)**(-4/3.), 'k--')

subplot(3,1,2, sharex=a)
for r in rates:
    plot(r[:,0], r[:,1].sum()-np.cumsum(r[:,1]))

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-1/3.), 'k--')
#plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
#plot(1e-2/np.log(np.arange(1000)), 'k:')

#ylim(0,0.025)
yscale('log'); #ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(3,1,3, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
a = subplot(3,1,1)
xscale('log'); yscale('log'); ylim(1e-3,0.025)
for r in smrates:
    plot(r[:,0], r[:,1])
plot(0.3*np.arange(1000)**(-4/3.), 'k--')

subplot(3,1,2, sharex=a)
for r in rates:
    plot(r[:,0], 1-np.cumsum(r[:,1]))

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-1/3.), 'k--')
#plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
#plot(1e-2/np.log(np.arange(1000)), 'k:')

#ylim(0,0.025)
yscale('log'); #ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(3,1,3, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
a = subplot(3,1,1)
xscale('log'); yscale('log'); ylim(1e-3,0.025)
for r in smrates:
    plot(r[:,0], r[:,1])
plot(0.3*np.arange(1000)**(-4/3.), 'k--')
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

subplot(3,1,2, sharex=a)
for r in rates:
    plot(r[:,0], 1-np.cumsum(r[:,1]))

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(10*np.arange(1000)**(-1/3.), 'k--')
#plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
#plot(1e-2/np.log(np.arange(1000)), 'k:')

#ylim(0,0.025)
yscale('log'); #ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(3,1,3, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
a = subplot(3,1,1)
xscale('log'); yscale('log'); ylim(1e-3,0.025)
for r in smrates:
    plot(r[:,0], r[:,1])
plot(0.3*np.arange(1000)**(-4/3.), 'k--')
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

subplot(3,1,2, sharex=a)
for r in rates:
    plot(r[:,0], 1-np.cumsum(r[:,1]))

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(np.arange(1000)**(-1/3.), 'k--')
#plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
#plot(1e-2/np.log(np.arange(1000)), 'k:')

#ylim(0,0.025)
yscale('log'); #ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(3,1,3, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
a = subplot(3,1,1)
xscale('log'); yscale('log'); ylim(1e-3,0.025)
for r in smrates:
    plot(r[:,0], r[:,1])
plot(0.3*np.arange(1000)**(-4/3.), 'k--')
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

subplot(3,1,2, sharex=a)
for r in rates:
    plot(r[:,0], 1-np.cumsum(r[:,1]))

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(3*np.arange(1000)**(-1/3.), 'k--')
#plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
#plot(1e-2/np.log(np.arange(1000)), 'k:')

#ylim(0,0.025)
yscale('log'); #ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(3,1,3, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
[l.get_color() for l in gca().lines]
for r in rates:
    #10 points per decade
    li = np.unique(np.logspace(1, np.log10(len(li))+1, 10*np.log10(len(li))))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
yscale('log')
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
for r in rates:
    #10 points per decade
    li = np.unique(np.logspace(1, np.log10(len(r))+1, 10*np.log10(len(r))))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
yscale('log')
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
for r in rates:
    #10 points per decade
    li = np.unique(np.logspace(1, np.log10(len(r))+1, 10*np.log10(len(r))).astype(int))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
yscale('log')
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
for r in rates:
    #10 points per decade
    li = np.unique(np.logspace(1, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
yscale('log')
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
for r in rates:
    #10 points per decade
    li = np.unique(np.logspace(1, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
yscale('log')
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
for r in rates:
    #10 points per decade
    li = np.unique(np.logspace(1, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
#yscale('log')
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
for r in rates:
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
#yscale('log')
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
for r in rates:
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
yscale('log')
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
for r in rates:
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])
plot(0.3*np.arange(1000)**(-4/3.), 'k--')
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
yscale('log')
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
for r in rates:
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])
plot(0.3*np.arange(1000)**(-4/3.), 'k--')
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
yscale('log')
ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
for r in rates:
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])
plot(0.3*np.arange(1000)**(-4/3.), 'k--')
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
yscale('log')
ylim(1e-4,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
for r in rates:
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])
plot(0.3*np.arange(1000)**(-4/3.), 'k--')
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
yscale('log')
ylim(3e-4,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
pa = '/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data'
ns = [
    '3_DenseGel/163A/1_percolation/1340_percolation',
    '2_MidGel/162B/1_percolation/1745',
    '1_DiluteGel/172A/1_percolation/172A_1206_percolation',
    #'1_DiluteGel/155C/1_percolation/155C_percolation_1645',
    '2_MidGel/170A/1_percolation/170A_1451_percolation',
    '2_MidGel/153A/1_percolation/153A_percolation_1540',
    '3_DenseGel/150A/1_percolation/150A_percolation_1227',
    '3_DenseGel/150A/1_percolation/150A_percolation_1239',
    '2_MidGel/162B/2_ageing/1815_ageing'
    '1_DiluteGel/172A/2_ageing/172A_1318_ageing'
    ]
xs = [
    xp.Experiment(os.path.join(pa, n+'.traj'))
    for n in ns
    ]
pa = '/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data'
ns = [
    '3_DenseGel/163A/1_percolation/1340_percolation',
    '2_MidGel/162B/1_percolation/1745',
    '1_DiluteGel/172A/1_percolation/172A_1206_percolation',
    #'1_DiluteGel/155C/1_percolation/155C_percolation_1645',
    '2_MidGel/170A/1_percolation/170A_1451_percolation',
    '2_MidGel/153A/1_percolation/153A_percolation_1540',
    '3_DenseGel/150A/1_percolation/150A_percolation_1227',
    '3_DenseGel/150A/1_percolation/150A_percolation_1239',
    '2_MidGel/162B/2_ageing/1815_ageing',
    '1_DiluteGel/172A/2_ageing/172A_1318_ageing'
    ]
xs = [
    xp.Experiment(os.path.join(pa, n+'.traj'))
    for n in ns
    ]
hists =[]
for x in xs:
    lengths = [
        np.loadtxt(name, int) 
        for t, name in x.enum('_broken', 'length') 
        if t>0]
    maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
    hists.append(np.array([
        np.histogram(l, np.arange(-1, maxlength+1))[0]
        for l in lengths]))
len(hists[2])
t0s = [68,151,33,64,200,104,0,0]
rates = [
    np.column_stack((
        np.arange(len(his)-t0),
        his[t0:,4:].sum(-1)*1./x.get_nb()[t0:-1]
        )) 
    for his, t0, x in zip(hists, t0s, xs)
    ]

#172A
rates[1] = np.vstack((rates[2], rates[-1]*[3,1]+[72*6-419,0]))
del rates[-1]
#162B
rates[1] = np.vstack((rates[1], rates[-1]*[1,1]+[35*6-159,0]))
del rates[-1]
#150A
rates[5] = np.vstack((rates[5], rates[6]+[8,0]))
del rates[6]
len(hists)
t0s = [68,151,33,64,200,104,0,0]
rates = [
    np.column_stack((
        np.arange(len(his)-t0),
        his[t0:,4:].sum(-1)*1./x.get_nb()[t0:-1]
        )) 
    for his, t0, x in zip(hists, t0s, xs)
    ]

#172A
rates[1] = np.vstack((rates[2], rates[-1]*[3,1]+[72*6-419,0]))
del rates[-1]
#162B
rates[1] = np.vstack((rates[1], rates[-1]*[1,1]+[35*6-159,0]))
del rates[-1]
#150A
rates[-2] = np.vstack((rates[-2], rates[-1]+[8,0]))
del rates[-1]
for r in rates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
#yscale('log')
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle.pdf')
t0s = [68,151,33,64,200,104,0,0]
rates = [
    np.column_stack((
        np.arange(len(his)-t0),
        his[t0:,4:].sum(-1)*1./x.get_nb()[t0:-1]
        )) 
    for his, t0, x in zip(hists, t0s, xs)
    ]

#172A
rates[2] = np.vstack((rates[2], rates[-1]*[3,1]+[72*6-419,0]))
del rates[-1]
#162B
rates[1] = np.vstack((rates[1], rates[-1]*[1,1]+[35*6-159,0]))
del rates[-1]
#150A
rates[-2] = np.vstack((rates[-2], rates[-1]+[8,0]))
del rates[-1]
for r in rates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
#yscale('log')
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle.pdf')
t0s = [68,151,33,64,200,104,0,0]
rates = [
    np.column_stack((
        np.arange(len(his)-t0),
        his[t0:,4:].sum(-1)*1./x.get_nb()[t0:-1]
        )) 
    for his, t0, x in zip(hists, t0s, xs)
    ]
print len(rates)
#172A
rates[2] = np.vstack((rates[2], rates[-1]*[3,1]+[72*6-419,0]))
del rates[-1]
print len(rates)
#162B
rates[1] = np.vstack((rates[1], rates[-1]*[1,1]+[35*6-159,0]))
del rates[-1]
print len(rates)
#150A
rates[-2] = np.vstack((rates[-2], rates[-1]+[8,0]))
del rates[-1]
print len(rates)
t0s = [68,151,33,64,200,104,0,0,0]
rates = [
    np.column_stack((
        np.arange(len(his)-t0),
        his[t0:,4:].sum(-1)*1./x.get_nb()[t0:-1]
        )) 
    for his, t0, x in zip(hists, t0s, xs)
    ]
print len(rates)
#172A
rates[2] = np.vstack((rates[2], rates[-1]*[3,1]+[72*6-419,0]))
del rates[-1]
print len(rates)
#162B
rates[1] = np.vstack((rates[1], rates[-1]*[1,1]+[35*6-159,0]))
del rates[-1]
print len(rates)
#150A
rates[-2] = np.vstack((rates[-2], rates[-1]+[8,0]))
del rates[-1]
print len(rates)
for r in rates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
#yscale('log')
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle.pdf')
t0 = 33
his = hists[2]
plot(np.arange(len(his)-t0), his)
his = hists[-1]
plot(np.arange(len(his)), his)
t0 = 33
his = hists[2]
plot(np.arange(len(his)-t0), his[t0:],4:].sum(-1))
his = hists[-1]
plot(np.arange(len(his)), his)
t0 = 33
his = hists[2]
plot(np.arange(len(his)-t0), his[t0:,4:].sum(-1))
his = hists[-1]
plot(np.arange(len(his)), his[:,4:].sum(-1))
t0 = 33
his = hists[2]
plot(np.arange(len(his)-t0), his[t0:,4:].sum(-1))
his = hists[-1]
plot(np.arange(len(his))+[72*6-419,0], his[:,4:].sum(-1))
t0 = 33
his = hists[2]
plot(np.arange(len(his)-t0), his[t0:,4:].sum(-1))
his = hists[-1]
plot(np.arange(len(his))+[72*6-419], his[:,4:].sum(-1))
t0 = 33
his = hists[2]
plot(np.arange(len(his)-t0), his[t0:,4:].sum(-1))
his = hists[-1]
plot(np.arange(len(his))+[419-72*6], his[:,4:].sum(-1))
t0 = 33
his = hists[2]
plot(np.arange(len(his)-t0), his[t0:,4:].sum(-1))
his = hists[-1]
plot(np.arange(len(his))+419, his[:,4:].sum(-1))
t0 = 33
his = hists[2]
plot(np.arange(len(his)-t0), his[t0:,4:].sum(-1))
his = hists[-1]
plot(np.arange(len(his))+419-33, his[:,4:].sum(-1))
t0 = 33
his = hists[2]
plot(np.arange(len(his)-t0), his[t0:,4:].sum(-1))
his = hists[-1]
plot(np.arange(len(his))+419-33, his[:,4:].sum(-1))
xlim(350,400)
t0 = 33
his = hists[2]
plot(np.arange(len(his)-t0), his[t0:,4:].sum(-1))
his = hists[-1]
plot(np.arange(len(his))+419-72+33, his[:,4:].sum(-1))
xlim(350,400)
t0 = 33
his = hists[2]
plot(np.arange(len(his)-t0), his[t0:,4:].sum(-1))
his = hists[-1]
plot(np.arange(len(his))+419-33, his[:,4:].sum(-1))
xlim(350,400)
t0s = [68,151,33,64,200,104,0,0,0]
rates = [
    np.column_stack((
        np.arange(len(his)-t0),
        his[t0:,4:].sum(-1)*1./x.get_nb()[t0:-1]
        )) 
    for his, t0, x in zip(hists, t0s, xs)
    ]
print len(rates)
#172A
rates[2] = np.vstack((rates[2], rates[-1]*[3,1]+[419-33,0]))
del rates[-1]
print len(rates)
#162B
rates[1] = np.vstack((rates[1], rates[-1]*[1,1]+[35*6-159,0]))
del rates[-1]
print len(rates)
#150A
rates[-2] = np.vstack((rates[-2], rates[-1]+[8,0]))
del rates[-1]
print len(rates)
for r in rates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
#yscale('log')
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle.pdf')
t0 = 33
his = hists[2]
plot(np.arange(len(his)-t0), his[t0:,4:].sum(-1))
his = hists[-1]
plot(np.arange(len(his))+419-33, his[:,4:].sum(-1))
#xlim(350,400)
t0 = 33
his = hists[2]
plot(np.arange(len(his)-t0), his[t0:,4:].sum(-1))
his = hists[-1][:400]
plot(np.arange(len(his))+419-33, his[:,4:].sum(-1))
#xlim(350,400)
t0 = 33
his = hists[2]
plot(np.arange(len(his)-t0), his[t0:,4:].sum(-1))
his = hists[-1][:450]
plot(np.arange(len(his))+419-33, his[:,4:].sum(-1))
#xlim(350,400)
t0 = 33
his = hists[2]
plot(np.arange(len(his)-t0), his[t0:,4:].sum(-1))
his = hists[-1][:500]
plot(np.arange(len(his))+419-33, his[:,4:].sum(-1))
#xlim(350,400)
t0 = 33
his = hists[2]
plot(np.arange(len(his)-t0), his[t0:,4:].sum(-1))
his = hists[-1][:480]
plot(np.arange(len(his))+419-33, his[:,4:].sum(-1))
#xlim(350,400)
t0 = 33
his = hists[2]
plot(np.arange(len(his)-t0), his[t0:,4:].sum(-1))
his = hists[-1][:460]
plot(np.arange(len(his))+419-33, his[:,4:].sum(-1))
#xlim(350,400)
t0 = 33
his = hists[2]
plot(np.arange(len(his)-t0), his[t0:,4:].sum(-1))
his = hists[-1][:470]
plot(np.arange(len(his))+419-33, his[:,4:].sum(-1))
#xlim(350,400)
t0 = 33
his = hists[2]
plot(np.arange(len(his)-t0), his[t0:,4:].sum(-1))
his = hists[-1][:465]
plot(np.arange(len(his))+419-33, his[:,4:].sum(-1))
#xlim(350,400)
t0 = 33
his = hists[2]
plot(np.arange(len(his)-t0), his[t0:,4:].sum(-1))
his = hists[-1][:460]
plot(np.arange(len(his))+419-33, his[:,4:].sum(-1))
#xlim(350,400)
t0s = [68,151,33,64,200,104,0,0,0]
rates = [
    np.column_stack((
        np.arange(len(his)-t0),
        his[t0:,4:].sum(-1)*1./x.get_nb()[t0:-1]
        )) 
    for his, t0, x in zip(hists, t0s, xs)
    ]
print len(rates)
#172A
rates[2] = np.vstack((rates[2], rates[-1][:460]*[3,1]+[419-33,0]))
del rates[-1]
print len(rates)
#162B
rates[1] = np.vstack((rates[1], rates[-1]*[1,1]+[35*6-159,0]))
del rates[-1]
print len(rates)
#150A
rates[-2] = np.vstack((rates[-2], rates[-1]+[8,0]))
del rates[-1]
print len(rates)
for r in rates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
#yscale('log')
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle.pdf')
get_ipython().system(u'ls -F --color /media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors/172A*_structure_factor.npy')
for r in rates:
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])
plot(0.3*np.arange(1000)**(-4/3.), 'k--')
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
yscale('log')
ylim(3e-4,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
qmaxs = []
names = [
    '163A_1340_percolation',
    '162B_percolation',
    '172A_1206_percolation',
    '170A_1451_percolation',
    '153A_percolation_1540',
    '150A_percolation_1227',
    '150A_percolation_1239',
    '162B_1815_ageing',
    '172A_1318_ageing'
    ]
imax = 44
for name, t0, x in zip(names, t0s, xs):
    Ss = np.load(os.path.join(
        '/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors',
        name+'_structure_factor.npy'
        ))
    q = np.arange(Ss.shape[1])*Ss[0,0]*2*x.radius
    qmaxs.append(
        np.column_stack((
            np.arange(len(Ss)-t0),
            (Ss[t0:,4:imax]*q[4:imax]).sum(1)/Ss[t0:,4:imax].sum(1)
            ))
        )
#172A
qmaxs[2] = np.vstack((qmaxs[2], qmaxs[-1][:460]*[3,1]+[419-33,0]))
del qmaxs[-1]
#162B
qmaxs[1] = np.vstack((qmaxs[1], qmaxs[-1]*[1,1]+[35*6-159,0]))
del qmaxs[-1]
#150A
qmaxs[5] = np.vstack((qmaxs[5], qmaxs[6]+[8,0]))
del qmaxs[6]
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
#plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
xscale('log');
yscale('log');
xlabel(r'$t/\tau_B$')
ylabel(r'$q_\mathrm{max}\sigma$')
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in rates:
    plot(r[:,0], r[:,1])
    
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
plot(0.3*np.arange(1000)**(-4/3.), 'k--')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
def smh(h, thr=100):
    """Smooth time dependant histograms by accumulation of at least thr counts. Returns times and counts."""
    dt=0
    ac = 0
    ts = []
    h2 = []
    for i, c in enumerate(h):
        if c > thr:
            if dt>0:
                h2.append(ac*1./dt)
                ts.append(i-dt) #note the starting time of the interval
                ac = 0
                dt = 0
            h2.append(c)
            ts.append(i)
        else:
            ac += c
            dt += 1
            if ac > thr or i+1==len(h):
                h2.append(ac*1./dt)
                ts.append(i+1-dt) #note the starting time of the interval
                ac = 0
                dt = 0
    return np.array(ts), np.array(h2)

smrates = [
    np.column_stack(smh(his[t0:,4:].sum(-1)))
    for his, t0 in zip(hists, t0s)
    ]
for r, x in zip(smrates, xs):
    r[:,1] /= x.get_nb()[t0:][r[:,0].astype(int)]

#172A
smrates[2] = np.vstack((smrates[2], smrates[-1][:460]*[3,1]+[419-33,0]))
del smrates[-1]
#162B
smrates[1] = np.vstack((smrates[1], smrates[-1]*[1,1]+[35*6-159,0]))
del smrates[-1]
#150A
smrates[5] = np.vstack((smrates[5], smrates[6]+[8,0]))
del smrates[6]
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
[l.get_color() for l in gca().lines]
def smh(h, thr=100):
    """Smooth time dependant histograms by accumulation of at least thr counts. Returns times and counts."""
    dt=0
    ac = 0
    ts = []
    h2 = []
    for i, c in enumerate(h):
        if c > thr:
            if dt>0:
                h2.append(ac*1./dt)
                ts.append(i-dt) #note the starting time of the interval
                ac = 0
                dt = 0
            h2.append(c)
            ts.append(i)
        else:
            ac += c
            dt += 1
            if ac > thr or i+1==len(h):
                h2.append(ac*1./dt)
                ts.append(i+1-dt) #note the starting time of the interval
                ac = 0
                dt = 0
    return np.array(ts), np.array(h2)

smrates = [
    np.column_stack(smh(his[t0:,4:].sum(-1)))
    for his, t0 in zip(hists, t0s)
    ]
for r, x in zip(smrates, xs):
    r[:,1] /= x.get_nb()[t0:][r[:,0].astype(int)]

#172A
smrates[2] = np.vstack((smrates[2], smrates[-1][:450]*[3,1]+[419-33,0]))
del smrates[-1]
#162B
smrates[1] = np.vstack((smrates[1], smrates[-1]*[1,1]+[35*6-159,0]))
del smrates[-1]
#150A
smrates[5] = np.vstack((smrates[5], smrates[6]+[8,0]))
del smrates[6]
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
[l.get_color() for l in gca().lines]
def smh(h, thr=100):
    """Smooth time dependant histograms by accumulation of at least thr counts. Returns times and counts."""
    dt=0
    ac = 0
    ts = []
    h2 = []
    for i, c in enumerate(h):
        if c > thr:
            if dt>0:
                h2.append(ac*1./dt)
                ts.append(i-dt) #note the starting time of the interval
                ac = 0
                dt = 0
            h2.append(c)
            ts.append(i)
        else:
            ac += c
            dt += 1
            if ac > thr or i+1==len(h):
                h2.append(ac*1./dt)
                ts.append(i+1-dt) #note the starting time of the interval
                ac = 0
                dt = 0
    return np.array(ts), np.array(h2)

smrates = [
    np.column_stack(smh(his[t0:,4:].sum(-1)))
    for his, t0 in zip(hists, t0s)
    ]
for r, x in zip(smrates, xs):
    r[:,1] /= x.get_nb()[t0:][r[:,0].astype(int)]

print len(smrates[-1])

#172A
smrates[2] = np.vstack((smrates[2], smrates[-1][:430]*[3,1]+[419-33,0]))
del smrates[-1]
#162B
smrates[1] = np.vstack((smrates[1], smrates[-1]*[1,1]+[35*6-159,0]))
del smrates[-1]
#150A
smrates[5] = np.vstack((smrates[5], smrates[6]+[8,0]))
del smrates[6]
def smh(h, thr=100):
    """Smooth time dependant histograms by accumulation of at least thr counts. Returns times and counts."""
    dt=0
    ac = 0
    ts = []
    h2 = []
    for i, c in enumerate(h):
        if c > thr:
            if dt>0:
                h2.append(ac*1./dt)
                ts.append(i-dt) #note the starting time of the interval
                ac = 0
                dt = 0
            h2.append(c)
            ts.append(i)
        else:
            ac += c
            dt += 1
            if ac > thr or i+1==len(h):
                h2.append(ac*1./dt)
                ts.append(i+1-dt) #note the starting time of the interval
                ac = 0
                dt = 0
    return np.array(ts), np.array(h2)

smrates = [
    np.column_stack(smh(his[t0:,4:].sum(-1)))
    for his, t0 in zip(hists, t0s)
    ]
for r, x in zip(smrates, xs):
    r[:,1] /= x.get_nb()[t0:][r[:,0].astype(int)]

print len(smrates[-1])

#172A
smrates[2] = np.vstack((smrates[2], smrates[-1][:-1]*[3,1]+[419-33,0]))
del smrates[-1]
#162B
smrates[1] = np.vstack((smrates[1], smrates[-1]*[1,1]+[35*6-159,0]))
del smrates[-1]
#150A
smrates[5] = np.vstack((smrates[5], smrates[6]+[8,0]))
del smrates[6]
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
[l.get_color() for l in gca().lines]
def smh(h, thr=100):
    """Smooth time dependant histograms by accumulation of at least thr counts. Returns times and counts."""
    dt=0
    ac = 0
    ts = []
    h2 = []
    for i, c in enumerate(h):
        if c > thr:
            if dt>0:
                h2.append(ac*1./dt)
                ts.append(i-dt) #note the starting time of the interval
                ac = 0
                dt = 0
            h2.append(c)
            ts.append(i)
        else:
            ac += c
            dt += 1
            if ac > thr or i+1==len(h):
                h2.append(ac*1./dt)
                ts.append(i+1-dt) #note the starting time of the interval
                ac = 0
                dt = 0
    return np.array(ts), np.array(h2)

smrates = [
    np.column_stack(smh(his[t0:,4:].sum(-1)))
    for his, t0 in zip(hists, t0s)
    ]
for r, x in zip(smrates, xs):
    r[:,1] /= x.get_nb()[t0:][r[:,0].astype(int)]

print len(smrates[-1])

#172A
smrates[2] = np.vstack((smrates[2], smrates[-1][:-2]*[3,1]+[419-33,0]))
del smrates[-1]
#162B
smrates[1] = np.vstack((smrates[1], smrates[-1]*[1,1]+[35*6-159,0]))
del smrates[-1]
#150A
smrates[5] = np.vstack((smrates[5], smrates[6]+[8,0]))
del smrates[6]
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
[l.get_color() for l in gca().lines]
def smh(h, thr=100):
    """Smooth time dependant histograms by accumulation of at least thr counts. Returns times and counts."""
    dt=0
    ac = 0
    ts = []
    h2 = []
    for i, c in enumerate(h):
        if c > thr:
            if dt>0:
                h2.append(ac*1./dt)
                ts.append(i-dt) #note the starting time of the interval
                ac = 0
                dt = 0
            h2.append(c)
            ts.append(i)
        else:
            ac += c
            dt += 1
            if ac > thr or i+1==len(h):
                h2.append(ac*1./dt)
                ts.append(i+1-dt) #note the starting time of the interval
                ac = 0
                dt = 0
    return np.array(ts), np.array(h2)

smrates = [
    np.column_stack(smh(his[t0:,4:].sum(-1)))
    for his, t0 in zip(hists, t0s)
    ]
for r, x in zip(smrates, xs):
    r[:,1] /= x.get_nb()[t0:][r[:,0].astype(int)]

print len(smrates[-1])

#172A
smrates[2] = np.vstack((smrates[2], smrates[-1][:-5]*[3,1]+[419-33,0]))
del smrates[-1]
#162B
smrates[1] = np.vstack((smrates[1], smrates[-1]*[1,1]+[35*6-159,0]))
del smrates[-1]
#150A
smrates[5] = np.vstack((smrates[5], smrates[6]+[8,0]))
del smrates[6]
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
[l.get_color() for l in gca().lines]
def smh(h, thr=100):
    """Smooth time dependant histograms by accumulation of at least thr counts. Returns times and counts."""
    dt=0
    ac = 0
    ts = []
    h2 = []
    for i, c in enumerate(h):
        if c > thr:
            if dt>0:
                h2.append(ac*1./dt)
                ts.append(i-dt) #note the starting time of the interval
                ac = 0
                dt = 0
            h2.append(c)
            ts.append(i)
        else:
            ac += c
            dt += 1
            if ac > thr or i+1==len(h):
                h2.append(ac*1./dt)
                ts.append(i+1-dt) #note the starting time of the interval
                ac = 0
                dt = 0
    return np.array(ts), np.array(h2)

smrates = [
    np.column_stack(smh(his[t0:,4:].sum(-1)))
    for his, t0 in zip(hists, t0s)
    ]
for r, x in zip(smrates, xs):
    r[:,1] /= x.get_nb()[t0:][r[:,0].astype(int)]

print len(smrates[-1])

#172A
smrates[2] = np.vstack((smrates[2], smrates[-1][:-4]*[3,1]+[419-33,0]))
del smrates[-1]
#162B
smrates[1] = np.vstack((smrates[1], smrates[-1]*[1,1]+[35*6-159,0]))
del smrates[-1]
#150A
smrates[5] = np.vstack((smrates[5], smrates[6]+[8,0]))
del smrates[6]
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
[l.get_color() for l in gca().lines]
def smh(h, thr=100):
    """Smooth time dependant histograms by accumulation of at least thr counts. Returns times and counts."""
    dt=0
    ac = 0
    ts = []
    h2 = []
    for i, c in enumerate(h):
        if c > thr:
            if dt>0:
                h2.append(ac*1./dt)
                ts.append(i-dt) #note the starting time of the interval
                ac = 0
                dt = 0
            h2.append(c)
            ts.append(i)
        else:
            ac += c
            dt += 1
            if ac > thr or i+1==len(h):
                h2.append(ac*1./dt)
                ts.append(i+1-dt) #note the starting time of the interval
                ac = 0
                dt = 0
    return np.array(ts), np.array(h2)

smrates = [
    np.column_stack(smh(his[t0:,4:].sum(-1)))
    for his, t0 in zip(hists, t0s)
    ]
for r, x in zip(smrates, xs):
    r[:,1] /= x.get_nb()[t0:][r[:,0].astype(int)]

print len(smrates[-1])

#172A
smrates[2] = np.vstack((smrates[2], smrates[-1][:-3]*[3,1]+[419-33,0]))
del smrates[-1]
#162B
smrates[1] = np.vstack((smrates[1], smrates[-1]*[1,1]+[35*6-159,0]))
del smrates[-1]
#150A
smrates[5] = np.vstack((smrates[5], smrates[6]+[8,0]))
del smrates[6]
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
[l.get_color() for l in gca().lines]
def smh(h, thr=100):
    """Smooth time dependant histograms by accumulation of at least thr counts. Returns times and counts."""
    dt=0
    ac = 0
    ts = []
    h2 = []
    for i, c in enumerate(h):
        if c > thr:
            if dt>0:
                h2.append(ac*1./dt)
                ts.append(i-dt) #note the starting time of the interval
                ac = 0
                dt = 0
            h2.append(c)
            ts.append(i)
        else:
            ac += c
            dt += 1
            if ac > thr or i+1==len(h):
                h2.append(ac*1./dt)
                ts.append(i+1-dt) #note the starting time of the interval
                ac = 0
                dt = 0
    return np.array(ts), np.array(h2)

smrates = [
    np.column_stack(smh(his[t0:,4:].sum(-1)))
    for his, t0 in zip(hists, t0s)
    ]
for r, x in zip(smrates, xs):
    r[:,1] /= x.get_nb()[t0:][r[:,0].astype(int)]

print len(smrates[-1])

#172A
smrates[2] = np.vstack((smrates[2], smrates[-1][:-4]*[3,1]+[419-33,0]))
del smrates[-1]
#162B
smrates[1] = np.vstack((smrates[1], smrates[-1]*[1,1]+[35*6-159,0]))
del smrates[-1]
#150A
smrates[5] = np.vstack((smrates[5], smrates[6]+[8,0]))
del smrates[6]
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
[l.get_color() for l in gca().lines]
figure(figsize=(6,12))
a = subplot(3,1,1)
xscale('log'); yscale('log'); ylim(1e-3,0.025)
for r in smrates:
    plot(r[:,0], r[:,1])
plot(0.3*np.arange(1000)**(-4/3.), 'k--')
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

subplot(3,1,2, sharex=a)
for r in rates:
    plot(r[:,0], 1-np.cumsum(r[:,1]))

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(3*np.arange(1000)**(-1/3.), 'k--')
#plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
#plot(1e-2/np.log(np.arange(1000)), 'k:')

#ylim(0,0.025)
yscale('log'); #ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(3,1,3, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
[l.get_color() for l in gca().lines]
pa = '/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data'
ns = [
    '3_DenseGel/163A/1_percolation/1340_percolation',
    '2_MidGel/162B/1_percolation/1745',
    '1_DiluteGel/172A/1_percolation/172A_1206_percolation',
    #'1_DiluteGel/155C/1_percolation/155C_percolation_1645',
    '2_MidGel/170A/1_percolation/170A_1451_percolation',
    '2_MidGel/153A/1_percolation/153A_percolation_1540',
    '3_DenseGel/150A/1_percolation/150A_percolation_1227',
    '3_DenseGel/150A/1_percolation/150A_percolation_1239',
    '2_MidGel/162B/2_ageing/1815_ageing',
    '1_DiluteGel/172A/2_ageing/172A_1318_ageing',
    '3_DenseGel/163A/2_ageing/1415_ageing',
    ]
xs = [
    xp.Experiment(os.path.join(pa, n+'.traj'))
    for n in ns
    ]
hists =[]
for x in xs:
    lengths = [
        np.loadtxt(name, int) 
        for t, name in x.enum('_broken', 'length') 
        if t>0]
    maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
    hists.append(np.array([
        np.histogram(l, np.arange(-1, maxlength+1))[0]
        for l in lengths]))
len(hists[-1])
len(hists[0])
t0 = 68
his = hists[0]
plot(np.arange(len(his)-t0), his[t0:,4:].sum(-1))
his = hists[-1][:460]
plot(np.arange(len(his))+205-68, his[:,4:].sum(-1))
#xlim(350,400)
t0 = 68
his = hists[0]
plot(np.arange(len(his)-t0), his[t0:,4:].sum(-1))
his = hists[-1][:460]
plot(np.arange(len(his))+205-68, his[:,4:].sum(-1))
xlim(130,150)
t0s = [68,151,33,64,200,104,0,0,0,0]
rates = [
    np.column_stack((
        np.arange(len(his)-t0),
        his[t0:,4:].sum(-1)*1./x.get_nb()[t0:-1]
        )) 
    for his, t0, x in zip(hists, t0s, xs)
    ]
print len(rates)
#163A
rates[0] = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
del rates[-1]
print len(rates)
#172A
rates[2] = np.vstack((rates[2], rates[-1][:460]*[3,1]+[419-33,0]))
del rates[-1]
print len(rates)
#162B
rates[1] = np.vstack((rates[1], rates[-1]*[1,1]+[35*6-159,0]))
del rates[-1]
print len(rates)
#150A
rates[-2] = np.vstack((rates[-2], rates[-1]+[8,0]))
del rates[-1]
print len(rates)
for r in rates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
#yscale('log')
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle.pdf')
for r in rates:
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])
plot(0.3*np.arange(1000)**(-4/3.), 'k--')
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
yscale('log')
ylim(3e-4,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
get_ipython().system(u'ls -F --color /media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors/163A*_structure_factor.npy')
get_ipython().system(u'ls -F --color /media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors/163A*_structure_factor.npy')
qmaxs = []
names = [
    '163A_1340_percolation',
    '162B_percolation',
    '172A_1206_percolation',
    '170A_1451_percolation',
    '153A_percolation_1540',
    '150A_percolation_1227',
    '150A_percolation_1239',
    '162B_1815_ageing',
    '172A_1318_ageing',
    '163A_1415_ageing',
    ]
imax = 44
for name, t0, x in zip(names, t0s, xs):
    Ss = np.load(os.path.join(
        '/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors',
        name+'_structure_factor.npy'
        ))
    q = np.arange(Ss.shape[1])*Ss[0,0]*2*x.radius
    qmaxs.append(
        np.column_stack((
            np.arange(len(Ss)-t0),
            (Ss[t0:,4:imax]*q[4:imax]).sum(1)/Ss[t0:,4:imax].sum(1)
            ))
        )
#172A
qmaxs[0] = np.vstack((qmaxs[0], qmaxs[-1][:460]*[3,1]+[205-68,0]))
del qmaxs[-1]
#172A
qmaxs[2] = np.vstack((qmaxs[2], qmaxs[-1][:460]*[3,1]+[419-33,0]))
del qmaxs[-1]
#162B
qmaxs[1] = np.vstack((qmaxs[1], qmaxs[-1]*[1,1]+[35*6-159,0]))
del qmaxs[-1]
#150A
qmaxs[5] = np.vstack((qmaxs[5], qmaxs[6]+[8,0]))
del qmaxs[6]
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
#plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
xscale('log');
yscale('log');
xlabel(r'$t/\tau_B$')
ylabel(r'$q_\mathrm{max}\sigma$')
qmaxs = []
names = [
    '163A_1340_percolation',
    '162B_percolation',
    '172A_1206_percolation',
    '170A_1451_percolation',
    '153A_percolation_1540',
    '150A_percolation_1227',
    '150A_percolation_1239',
    '162B_1815_ageing',
    '172A_1318_ageing',
    '163A_1415_ageing',
    ]
imax = 44
for name, t0, x in zip(names, t0s, xs):
    Ss = np.load(os.path.join(
        '/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors',
        name+'_structure_factor.npy'
        ))
    q = np.arange(Ss.shape[1])*Ss[0,0]*2*x.radius
    qmaxs.append(
        np.column_stack((
            np.arange(len(Ss)-t0),
            (Ss[t0:,4:imax]*q[4:imax]).sum(1)/Ss[t0:,4:imax].sum(1)
            ))
        )
#172A
qmaxs[0] = np.vstack((qmaxs[0], qmaxs[-1]*[3,1]+[205-68,0]))
del qmaxs[-1]
#172A
qmaxs[2] = np.vstack((qmaxs[2], qmaxs[-1][:460]*[3,1]+[419-33,0]))
del qmaxs[-1]
#162B
qmaxs[1] = np.vstack((qmaxs[1], qmaxs[-1]*[1,1]+[35*6-159,0]))
del qmaxs[-1]
#150A
qmaxs[5] = np.vstack((qmaxs[5], qmaxs[6]+[8,0]))
del qmaxs[6]
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
#plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
xscale('log');
yscale('log');
xlabel(r'$t/\tau_B$')
ylabel(r'$q_\mathrm{max}\sigma$')
qmaxs = []
names = [
    '163A_1340_percolation',
    '162B_percolation',
    '172A_1206_percolation',
    '170A_1451_percolation',
    '153A_percolation_1540',
    '150A_percolation_1227',
    '150A_percolation_1239',
    '162B_1815_ageing',
    '172A_1318_ageing',
    '163A_1415_ageing',
    ]
imax = 44
for name, t0, x in zip(names, t0s, xs):
    Ss = np.load(os.path.join(
        '/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors',
        name+'_structure_factor.npy'
        ))
    q = np.arange(Ss.shape[1])*Ss[0,0]*2*x.radius
    qmaxs.append(
        np.column_stack((
            np.arange(len(Ss)-t0),
            (Ss[t0:,4:imax]*q[4:imax]).sum(1)/Ss[t0:,4:imax].sum(1)
            ))
        )
#172A
qmaxs[0] = np.vstack((qmaxs[0], qmaxs[-1]*[1,1]+[205-68,0]))
del qmaxs[-1]
#172A
qmaxs[2] = np.vstack((qmaxs[2], qmaxs[-1][:460]*[3,1]+[419-33,0]))
del qmaxs[-1]
#162B
qmaxs[1] = np.vstack((qmaxs[1], qmaxs[-1]*[1,1]+[35*6-159,0]))
del qmaxs[-1]
#150A
qmaxs[5] = np.vstack((qmaxs[5], qmaxs[6]+[8,0]))
del qmaxs[6]
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
#plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
xscale('log');
yscale('log');
xlabel(r'$t/\tau_B$')
ylabel(r'$q_\mathrm{max}\sigma$')
qmaxs = []
names = [
    '163A_1340_percolation',
    '162B_percolation',
    '172A_1206_percolation',
    '170A_1451_percolation',
    '153A_percolation_1540',
    '150A_percolation_1227',
    '150A_percolation_1239',
    '162B_1815_ageing',
    '172A_1318_ageing',
    '163A_1415_ageing',
    ]
imax = 44
for name, t0, x in zip(names, t0s, xs):
    Ss = np.load(os.path.join(
        '/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors',
        name+'_structure_factor.npy'
        ))
    q = np.arange(Ss.shape[1])*Ss[0,0]*2*x.radius
    qmaxs.append(
        np.column_stack((
            np.arange(len(Ss)-t0),
            (Ss[t0:,4:imax]*q[4:imax]).sum(1)/Ss[t0:,4:imax].sum(1)
            ))
        )
#172A
qmaxs[0] = np.vstack((qmaxs[0], qmaxs[-1]*[3,1]+[205-68,0]))
del qmaxs[-1]
#172A
qmaxs[2] = np.vstack((qmaxs[2], qmaxs[-1][:460]*[3,1]+[419-33,0]))
del qmaxs[-1]
#162B
qmaxs[1] = np.vstack((qmaxs[1], qmaxs[-1]*[1,1]+[35*6-159,0]))
del qmaxs[-1]
#150A
qmaxs[5] = np.vstack((qmaxs[5], qmaxs[6]+[8,0]))
del qmaxs[6]
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
#plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
xscale('log');
yscale('log');
xlabel(r'$t/\tau_B$')
ylabel(r'$q_\mathrm{max}\sigma$')
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in rates:
    plot(r[:,0], r[:,1])
    
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
plot(0.3*np.arange(1000)**(-4/3.), 'k--')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
def smh(h, thr=100):
    """Smooth time dependant histograms by accumulation of at least thr counts. Returns times and counts."""
    dt=0
    ac = 0
    ts = []
    h2 = []
    for i, c in enumerate(h):
        if c > thr:
            if dt>0:
                h2.append(ac*1./dt)
                ts.append(i-dt) #note the starting time of the interval
                ac = 0
                dt = 0
            h2.append(c)
            ts.append(i)
        else:
            ac += c
            dt += 1
            if ac > thr or i+1==len(h):
                h2.append(ac*1./dt)
                ts.append(i+1-dt) #note the starting time of the interval
                ac = 0
                dt = 0
    return np.array(ts), np.array(h2)

smrates = [
    np.column_stack(smh(his[t0:,4:].sum(-1)))
    for his, t0 in zip(hists, t0s)
    ]
for r, x in zip(smrates, xs):
    r[:,1] /= x.get_nb()[t0:][r[:,0].astype(int)]

print len(smrates[-1])
#163A
smrates[0] = np.vstack((smrates[0], smrates[-1]*[3,1]+[205-68,0]))
del smrates[-1]
#172A
smrates[2] = np.vstack((smrates[2], smrates[-1][:-4]*[3,1]+[419-33,0]))
del smrates[-1]
#162B
smrates[1] = np.vstack((smrates[1], smrates[-1]*[1,1]+[35*6-159,0]))
del smrates[-1]
#150A
smrates[5] = np.vstack((smrates[5], smrates[6]+[8,0]))
del smrates[6]
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
[l.get_color() for l in gca().lines]
t0s = [68,151,33,64,200,104,0,0,0,0]
rates = [
    np.column_stack((
        np.arange(len(his)-t0),
        his[t0:,4:].sum(-1)*1./x.get_nb()[t0:-1]
        )) 
    for his, t0, x in zip(hists, t0s, xs)
    ]
print len(rates)
#163A
rates[0] = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
del rates[-1]
print len(rates)
#172A
rates[2] = np.vstack((rates[2], rates[-1][:460]*[3,1]+[419-33,0]))
del rates[-1]
print len(rates)
#162B
rates[1] = np.vstack((rates[1], rates[-1][:-1]*[1,1]+[35*6-159,0]))
del rates[-1]
print len(rates)
#150A
rates[-2] = np.vstack((rates[-2], rates[-1]+[8,0]))
del rates[-1]
print len(rates)
for r in rates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
#yscale('log')
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle.pdf')
t0s = [68,151,33,64,200,104,0,0,0,0]
rates = [
    np.column_stack((
        np.arange(len(his)-t0),
        his[t0:,4:].sum(-1)*1./x.get_nb()[t0:-1]
        )) 
    for his, t0, x in zip(hists, t0s, xs)
    ]
print len(rates)
#163A
rates[0] = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
del rates[-1]
print len(rates)
#172A
rates[2] = np.vstack((rates[2], rates[-1][:460]*[3,1]+[419-33,0]))
del rates[-1]
print len(rates)
#162B
rates[1] = np.vstack((rates[1], rates[-1][:-5]*[1,1]+[35*6-159,0]))
del rates[-1]
print len(rates)
#150A
rates[-2] = np.vstack((rates[-2], rates[-1]+[8,0]))
del rates[-1]
print len(rates)
for r in rates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
#yscale('log')
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle.pdf')
t0 = 151
his = hists[1]
plot(np.arange(len(his)-t0), his[t0:,4:].sum(-1))
his = hists[7]
plot(np.arange(len(his))+35*6-159, his[:,4:].sum(-1))
t0s = [68,151,33,64,200,104,0,0,0,0]
rates = [
    np.column_stack((
        np.arange(len(his)-t0),
        his[t0:,4:].sum(-1)*1./x.get_nb()[t0:-1]
        )) 
    for his, t0, x in zip(hists, t0s, xs)
    ]
print len(rates)
#163A
rates[0] = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
del rates[-1]
print len(rates)
#172A
rates[2] = np.vstack((rates[2], rates[-1][:460]*[3,1]+[419-33,0]))
del rates[-1]
print len(rates)
#162B
rates[1] = np.vstack((rates[1], rates[-1]*[1,1]+[35*6-159,0]))
del rates[-1]
print len(rates)
#150A
rates[-2] = np.vstack((rates[-2], rates[-1]+[8,0]))
del rates[-1]
print len(rates)
figure(figsize=(6,12))
a = subplot(3,1,1)
xscale('log'); yscale('log'); ylim(1e-3,0.025)
for r in smrates:
    plot(r[:,0], r[:,1])
plot(0.3*np.arange(1000)**(-4/3.), 'k--')
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

subplot(3,1,2, sharex=a)
for r in rates:
    plot(r[:,0], 1-np.cumsum(r[:,1]))

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(3*np.arange(1000)**(-1/3.), 'k--')
#plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
#plot(1e-2/np.log(np.arange(1000)), 'k:')

#ylim(0,0.025)
yscale('log'); #ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(3,1,3, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
[l.get_color() for l in gca().lines]
pa = '/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data'
ns = [
    '3_DenseGel/163A/1_percolation/1340_percolation',
    '2_MidGel/162B/1_percolation/1745',
    '1_DiluteGel/172A/1_percolation/172A_1206_percolation',
    #'1_DiluteGel/155C/1_percolation/155C_percolation_1645',
    '2_MidGel/170A/1_percolation/170A_1451_percolation',
    '2_MidGel/153A/1_percolation/153A_percolation_1540',
    '3_DenseGel/150A/1_percolation/150A_percolation_1227',
    '3_DenseGel/150A/1_percolation/150A_percolation_1239',
    '2_MidGel/162B/2_ageing/1815_ageing',
    '1_DiluteGel/172A/2_ageing/172A_1318_ageing',
    '3_DenseGel/163A/2_ageing/1415_ageing',
    '2_MidGel/153A/2_percolation/153A_percolation_1610',
    '2_MidGel/153A/3_ageing/153A_ageing_1630',
    ]
xs = [
    xp.Experiment(os.path.join(pa, n+'.traj'))
    for n in ns
    ]
hists =[]
for x in xs:
    lengths = [
        np.loadtxt(name, int) 
        for t, name in x.enum('_broken', 'length') 
        if t>0]
    maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
    hists.append(np.array([
        np.histogram(l, np.arange(-1, maxlength+1))[0]
        for l in lengths]))
len(hists[4])
t0 = 64
his = hists[4]
plot(np.arange(len(his)-t0), his[t0:,4:].sum(-1))
his = hists[10]
plot(np.arange(len(his))+64-300, his[:,4:].sum(-1))
t0 = 64
his = hists[4]
plot(np.arange(len(his)-t0), his[t0:,4:].sum(-1))
his = hists[10]
plot(np.arange(len(his))+300-64, his[:,4:].sum(-1))
len(hists[4])+len(hist[10])
len(hists[4])+len(hists[10])
t0 = 64
his = hists[4]
plot(np.arange(len(his)-t0), his[t0:,4:].sum(-1))
his = hists[10]
plot(np.arange(len(his))+300-64, his[:,4:].sum(-1))
his = hists[11]
plot(np.arange(len(his))+300+167-64, his[:,4:].sum(-1))
t0 = 64
his = hists[4]
plot(np.arange(len(his)-t0), his[t0:,4:].sum(-1))
his = hists[10]
plot(np.arange(len(his))+300-64, his[:,4:].sum(-1))
his = hists[11]
plot(np.arange(len(his))+300+167-64, his[:,4:].sum(-1)/3)
t0 = 64
his = hists[4]
plot(np.arange(len(his)-t0), his[t0:,4:].sum(-1))
his = hists[10]
plot(np.arange(len(his))+300-64, his[:,4:].sum(-1))
his = hists[11]
plot(np.arange(len(his))+300+167-64, his[:,4:].sum(-1)/4)
t0 = 64
his = hists[4]
plot(np.arange(len(his)-t0), his[t0:,4:].sum(-1))
his = hists[10]
plot(np.arange(len(his))+300-64, his[:,4:].sum(-1))
his = hists[11]
plot(np.arange(len(his))+300+167-64, his[:,4:].sum(-1)/12)
t0s = [68,151,33,64,200,104] + [0]*len(hists)
rates = [
    np.column_stack((
        np.arange(len(his)-t0),
        his[t0:,4:].sum(-1)*1./x.get_nb()[t0:-1]
        )) 
    for his, t0, x in zip(hists, t0s, xs)
    ]
print len(rates)
#153A
rates[4] = np.vstack((rates[0], rates[-2]*[1,1]+[300-64,0], rates[-1]*[3,1]+[467-64,0]))
del rates[-2:]
#163A
rates[0] = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
del rates[-1]
print len(rates)
#172A
rates[2] = np.vstack((rates[2], rates[-1][:460]*[3,1]+[419-33,0]))
del rates[-1]
print len(rates)
#162B
rates[1] = np.vstack((rates[1], rates[-1]*[1,1]+[35*6-159,0]))
del rates[-1]
print len(rates)
#150A
rates[-2] = np.vstack((rates[-2], rates[-1]+[8,0]))
del rates[-1]
print len(rates)
t0s = [68,151,33,64,200,104] + [0]*len(hists)
rates = [
    np.column_stack((
        np.arange(len(his)-t0),
        his[t0:,4:].sum(-1)*1./x.get_nb()[t0:-1]
        )) 
    for his, t0, x in zip(hists, t0s, xs)
    ]
print len(rates)
#153A
rates[4] = np.vstack((rates[0], rates[-2]*[1,1]+[300-64,0], rates[-1]*[3,1]+[467-64,0]))
del rates[-2:]
print len(rates)
#163A
rates[0] = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
del rates[-1]
print len(rates)
#172A
rates[2] = np.vstack((rates[2], rates[-1][:460]*[3,1]+[419-33,0]))
del rates[-1]
print len(rates)
#162B
rates[1] = np.vstack((rates[1], rates[-1]*[1,1]+[35*6-159,0]))
del rates[-1]
print len(rates)
#150A
rates[-2] = np.vstack((rates[-2], rates[-1]+[8,0]))
del rates[-1]
print len(rates)
for r in rates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
#yscale('log')
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle.pdf')
for r in rates:
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])
plot(0.3*np.arange(1000)**(-4/3.), 'k--')
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
yscale('log')
ylim(3e-4,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
get_ipython().system(u'ls -F --color /media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors/153A*_structure_factor.npy')
qmaxs = []
names = [
    '163A_1340_percolation',
    '162B_percolation',
    '172A_1206_percolation',
    '170A_1451_percolation',
    '153A_percolation_1540',
    '150A_percolation_1227',
    '150A_percolation_1239',
    '162B_1815_ageing',
    '172A_1318_ageing',
    '163A_1415_ageing',
    '153A_percolation_1610',
    '153A_ageing_1630',
    ]
imax = 44
for name, t0, x in zip(names, t0s, xs):
    Ss = np.load(os.path.join(
        '/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors',
        name+'_structure_factor.npy'
        ))
    q = np.arange(Ss.shape[1])*Ss[0,0]*2*x.radius
    qmaxs.append(
        np.column_stack((
            np.arange(len(Ss)-t0),
            (Ss[t0:,4:imax]*q[4:imax]).sum(1)/Ss[t0:,4:imax].sum(1)
            ))
        )
#153A
qmaxs[4] = np.vstack((qmaxs[0], qmaxs[-2]*[1,1]+[300-64,0], qmaxs[-1]*[3,1]+[467-64,0]))
del qmaxs[-2:]
#172A
qmaxs[0] = np.vstack((qmaxs[0], qmaxs[-1]*[3,1]+[205-68,0]))
del qmaxs[-1]
#172A
qmaxs[2] = np.vstack((qmaxs[2], qmaxs[-1][:460]*[3,1]+[419-33,0]))
del qmaxs[-1]
#162B
qmaxs[1] = np.vstack((qmaxs[1], qmaxs[-1]*[1,1]+[35*6-159,0]))
del qmaxs[-1]
#150A
qmaxs[5] = np.vstack((qmaxs[5], qmaxs[6]+[8,0]))
del qmaxs[6]
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
#plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
xscale('log');
yscale('log');
xlabel(r'$t/\tau_B$')
ylabel(r'$q_\mathrm{max}\sigma$')
qmaxs = []
names = [
    '163A_1340_percolation',
    '162B_percolation',
    '172A_1206_percolation',
    '170A_1451_percolation',
    '153A_percolation_1540',
    '150A_percolation_1227',
    '150A_percolation_1239',
    '162B_1815_ageing',
    '172A_1318_ageing',
    '163A_1415_ageing',
    '153A_percolation_1610',
    '153A_ageing_1630',
    ]
imax = 44
for name, t0, x in zip(names, t0s, xs):
    Ss = np.load(os.path.join(
        '/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors',
        name+'_structure_factor.npy'
        ))
    q = np.arange(Ss.shape[1])*Ss[0,0]*2*x.radius
    qmaxs.append(
        np.column_stack((
            np.arange(len(Ss)-t0),
            (Ss[t0:,4:imax]*q[4:imax]).sum(1)/Ss[t0:,4:imax].sum(1)
            ))
        )
#153A
qmaxs[4] = np.vstack((qmaxs[0], qmaxs[-2]*[1,0.5]+[300-64,0], qmaxs[-1]*[3,1]+[467-64,0]))
del qmaxs[-2:]
#172A
qmaxs[0] = np.vstack((qmaxs[0], qmaxs[-1]*[3,1]+[205-68,0]))
del qmaxs[-1]
#172A
qmaxs[2] = np.vstack((qmaxs[2], qmaxs[-1][:460]*[3,1]+[419-33,0]))
del qmaxs[-1]
#162B
qmaxs[1] = np.vstack((qmaxs[1], qmaxs[-1]*[1,1]+[35*6-159,0]))
del qmaxs[-1]
#150A
qmaxs[5] = np.vstack((qmaxs[5], qmaxs[6]+[8,0]))
del qmaxs[6]
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
#plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
xscale('log');
yscale('log');
xlabel(r'$t/\tau_B$')
ylabel(r'$q_\mathrm{max}\sigma$')
qmaxs = []
names = [
    '163A_1340_percolation',
    '162B_percolation',
    '172A_1206_percolation',
    '170A_1451_percolation',
    '153A_percolation_1540',
    '150A_percolation_1227',
    '150A_percolation_1239',
    '162B_1815_ageing',
    '172A_1318_ageing',
    '163A_1415_ageing',
    '153A_percolation_1610',
    '153A_ageing_1630',
    ]
imax = 44
for name, t0, x in zip(names, t0s, xs):
    Ss = np.load(os.path.join(
        '/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors',
        name+'_structure_factor.npy'
        ))
    q = np.arange(Ss.shape[1])*Ss[0,0]*2*x.radius
    qmaxs.append(
        np.column_stack((
            np.arange(len(Ss)-t0),
            (Ss[t0:,4:imax]*q[4:imax]).sum(1)/Ss[t0:,4:imax].sum(1)
            ))
        )
#153A
qmaxs[4] = np.vstack((qmaxs[0], qmaxs[-2]*[1,0.25]+[300-64,0], qmaxs[-1]*[3,1]+[467-64,0]))
del qmaxs[-2:]
#172A
qmaxs[0] = np.vstack((qmaxs[0], qmaxs[-1]*[3,1]+[205-68,0]))
del qmaxs[-1]
#172A
qmaxs[2] = np.vstack((qmaxs[2], qmaxs[-1][:460]*[3,1]+[419-33,0]))
del qmaxs[-1]
#162B
qmaxs[1] = np.vstack((qmaxs[1], qmaxs[-1]*[1,1]+[35*6-159,0]))
del qmaxs[-1]
#150A
qmaxs[5] = np.vstack((qmaxs[5], qmaxs[6]+[8,0]))
del qmaxs[6]
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
#plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
xscale('log');
yscale('log');
xlabel(r'$t/\tau_B$')
ylabel(r'$q_\mathrm{max}\sigma$')
qmaxs = []
names = [
    '163A_1340_percolation',
    '162B_percolation',
    '172A_1206_percolation',
    '170A_1451_percolation',
    '153A_percolation_1540',
    '150A_percolation_1227',
    '150A_percolation_1239',
    '162B_1815_ageing',
    '172A_1318_ageing',
    '163A_1415_ageing',
    '153A_percolation_1610',
    '153A_ageing_1630',
    ]
imax = 44
for name, t0, x in zip(names, t0s, xs):
    Ss = np.load(os.path.join(
        '/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors',
        name+'_structure_factor.npy'
        ))
    q = np.arange(Ss.shape[1])*Ss[0,0]*2*x.radius
    qmaxs.append(
        np.column_stack((
            np.arange(len(Ss)-t0),
            (Ss[t0:,4:imax]*q[4:imax]).sum(1)/Ss[t0:,4:imax].sum(1)
            ))
        )
#153A
qmaxs[4] = np.vstack((qmaxs[0], qmaxs[-2]*[1,1/3.]+[300-64,0], qmaxs[-1]*[3,1]+[467-64,0]))
del qmaxs[-2:]
#172A
qmaxs[0] = np.vstack((qmaxs[0], qmaxs[-1]*[3,1]+[205-68,0]))
del qmaxs[-1]
#172A
qmaxs[2] = np.vstack((qmaxs[2], qmaxs[-1][:460]*[3,1]+[419-33,0]))
del qmaxs[-1]
#162B
qmaxs[1] = np.vstack((qmaxs[1], qmaxs[-1]*[1,1]+[35*6-159,0]))
del qmaxs[-1]
#150A
qmaxs[5] = np.vstack((qmaxs[5], qmaxs[6]+[8,0]))
del qmaxs[6]
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
#plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
xscale('log');
yscale('log');
xlabel(r'$t/\tau_B$')
ylabel(r'$q_\mathrm{max}\sigma$')
qmaxs = []
names = [
    '163A_1340_percolation',
    '162B_percolation',
    '172A_1206_percolation',
    '170A_1451_percolation',
    '153A_percolation_1540',
    '150A_percolation_1227',
    '150A_percolation_1239',
    '162B_1815_ageing',
    '172A_1318_ageing',
    '163A_1415_ageing',
    '153A_percolation_1610',
    '153A_ageing_1630',
    ]
imax = 44
for name, t0, x in zip(names, t0s, xs):
    Ss = np.load(os.path.join(
        '/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors',
        name+'_structure_factor.npy'
        ))
    q = np.arange(Ss.shape[1])*Ss[0,0]*2*x.radius
    qmaxs.append(
        np.column_stack((
            np.arange(len(Ss)-t0),
            (Ss[t0:,4:imax]*q[4:imax]).sum(1)/Ss[t0:,4:imax].sum(1)
            ))
        )
#153A
qmaxs[4] = np.vstack((qmaxs[0], qmaxs[-2]*[1,1/2.]+[300-64,0], qmaxs[-1]*[3,1]+[467-64,0]))
del qmaxs[-2:]
#172A
qmaxs[0] = np.vstack((qmaxs[0], qmaxs[-1]*[3,1]+[205-68,0]))
del qmaxs[-1]
#172A
qmaxs[2] = np.vstack((qmaxs[2], qmaxs[-1][:460]*[3,1]+[419-33,0]))
del qmaxs[-1]
#162B
qmaxs[1] = np.vstack((qmaxs[1], qmaxs[-1]*[1,1]+[35*6-159,0]))
del qmaxs[-1]
#150A
qmaxs[5] = np.vstack((qmaxs[5], qmaxs[6]+[8,0]))
del qmaxs[6]
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
#plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
xscale('log');
yscale('log');
xlabel(r'$t/\tau_B$')
ylabel(r'$q_\mathrm{max}\sigma$')
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in rates:
    plot(r[:,0], r[:,1])
    
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
plot(0.3*np.arange(1000)**(-4/3.), 'k--')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
def smh(h, thr=100):
    """Smooth time dependant histograms by accumulation of at least thr counts. Returns times and counts."""
    dt=0
    ac = 0
    ts = []
    h2 = []
    for i, c in enumerate(h):
        if c > thr:
            if dt>0:
                h2.append(ac*1./dt)
                ts.append(i-dt) #note the starting time of the interval
                ac = 0
                dt = 0
            h2.append(c)
            ts.append(i)
        else:
            ac += c
            dt += 1
            if ac > thr or i+1==len(h):
                h2.append(ac*1./dt)
                ts.append(i+1-dt) #note the starting time of the interval
                ac = 0
                dt = 0
    return np.array(ts), np.array(h2)

smrates = [
    np.column_stack(smh(his[t0:,4:].sum(-1)))
    for his, t0 in zip(hists, t0s)
    ]
for r, x in zip(smrates, xs):
    r[:,1] /= x.get_nb()[t0:][r[:,0].astype(int)]

print len(smrates[-1])
#153A
smrates[4] = np.vstack((smrates[0], smrates[-2]*[1,1]+[300-64,0], smrates[-1]*[3,1]+[467-64,0]))
del smrates[-2:]
#163A
smrates[0] = np.vstack((smrates[0], smrates[-1]*[3,1]+[205-68,0]))
del smrates[-1]
#172A
smrates[2] = np.vstack((smrates[2], smrates[-1][:-4]*[3,1]+[419-33,0]))
del smrates[-1]
#162B
smrates[1] = np.vstack((smrates[1], smrates[-1]*[1,1]+[35*6-159,0]))
del smrates[-1]
#150A
smrates[5] = np.vstack((smrates[5], smrates[6]+[8,0]))
del smrates[6]
figure(figsize=(6,12))
a = subplot(2,1,1)
for r in smrates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)

plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(5e-3*np.arange(1000)**(-1/4.), 'k-.')
plot(1e-2/np.log(np.arange(1000)), 'k:')
xscale('log')
#ylim(0,0.025)
yscale('log'); ylim(1e-3,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')

subplot(2,1,2, sharex=a)
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
yscale('log');
[l.get_color() for l in gca().lines]
for i,L in enumerate([3,4,5,6]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,4:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x in zip(hists, t0s, xs)
        ]

    #153A
    rates[4] = np.vstack((rates[0], rates[-2]*[1,1]+[300-64,0], rates[-1]*[3,1]+[467-64,0]))
    del rates[-2:]

    #163A
    rates[0] = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    del rates[-1]

    #172A
    rates[2] = np.vstack((rates[2], rates[-1][:460]*[3,1]+[419-33,0]))
    del rates[-1]

    #162B
    rates[1] = np.vstack((rates[1], rates[-1]*[1,1]+[35*6-159,0]))
    del rates[-1]

    #150A
    rates[-2] = np.vstack((rates[-2], rates[-1]+[8,0]))
    del rates[-1]
    
    subplot(2,2,i+1)
    for r in rates:
        #10 points per decade
        li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
        li = np.concatenate(([0],li))
        plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])
    plot(0.3*np.arange(1000)**(-4/3.), 'k--')
    for t, c in zip([5,13,148], 'bgr'): 
        axvline(t, color=c)
    xscale('log')
    yscale('log')
    ylim(3e-4,0.025)
    xlim(1,1e3)
    xlabel(r'$t/\tau_B$')
    ylabel(r'strand rupture per particle')
figure(figsize=(12,12))
for i,L in enumerate([3,4,5,6]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x in zip(hists, t0s, xs)
        ]

    #153A
    rates[4] = np.vstack((rates[0], rates[-2]*[1,1]+[300-64,0], rates[-1]*[3,1]+[467-64,0]))
    del rates[-2:]

    #163A
    rates[0] = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    del rates[-1]

    #172A
    rates[2] = np.vstack((rates[2], rates[-1][:460]*[3,1]+[419-33,0]))
    del rates[-1]

    #162B
    rates[1] = np.vstack((rates[1], rates[-1]*[1,1]+[35*6-159,0]))
    del rates[-1]

    #150A
    rates[-2] = np.vstack((rates[-2], rates[-1]+[8,0]))
    del rates[-1]
    
    subplot(2,2,i+1)
    for r in rates:
        #10 points per decade
        li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
        li = np.concatenate(([0],li))
        plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])
    plot(0.3*np.arange(1000)**(-4/3.), 'k--')
    for t, c in zip([5,13,148], 'bgr'): 
        axvline(t, color=c)
    xscale('log')
    yscale('log')
    ylim(3e-4,0.025)
    xlim(1,1e3)
    xlabel(r'$t/\tau_B$')
    ylabel(r'strand rupture per particle')
figure(figsize=(12,12))
for i,L in enumerate([3,4,5,6]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x in zip(hists, t0s, xs)
        ]

    #153A
    rates[4] = np.vstack((rates[0], rates[-2]*[1,1]+[300-64,0], rates[-1]*[3,1]+[467-64,0]))
    del rates[-2:]

    #163A
    rates[0] = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    del rates[-1]

    #172A
    rates[2] = np.vstack((rates[2], rates[-1][:460]*[3,1]+[419-33,0]))
    del rates[-1]

    #162B
    rates[1] = np.vstack((rates[1], rates[-1]*[1,1]+[35*6-159,0]))
    del rates[-1]

    #150A
    rates[-2] = np.vstack((rates[-2], rates[-1]+[8,0]))
    del rates[-1]
    
    subplot(2,2,i+1)
    for r in rates:
        #10 points per decade
        li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
        li = np.concatenate(([0],li))
        plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])
    plot(0.3*np.arange(1000)**(-4/3.), 'k--')
    for t, c in zip([5,13,148], 'bgr'): 
        axvline(t, color=c)
    xscale('log')
    yscale('log')
    ylim(1e-4,1e-1)
    xlim(1,1e3)
    xlabel(r'$t/\tau_B$')
    ylabel(r'strand rupture per particle')
figure(figsize=(12,12))
for i,L in enumerate([3,4,8,16]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x in zip(hists, t0s, xs)
        ]

    #153A
    rates[4] = np.vstack((rates[0], rates[-2]*[1,1]+[300-64,0], rates[-1]*[3,1]+[467-64,0]))
    del rates[-2:]

    #163A
    rates[0] = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    del rates[-1]

    #172A
    rates[2] = np.vstack((rates[2], rates[-1][:460]*[3,1]+[419-33,0]))
    del rates[-1]

    #162B
    rates[1] = np.vstack((rates[1], rates[-1]*[1,1]+[35*6-159,0]))
    del rates[-1]

    #150A
    rates[-2] = np.vstack((rates[-2], rates[-1]+[8,0]))
    del rates[-1]
    
    subplot(2,2,i+1)
    for r in rates:
        #10 points per decade
        li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
        li = np.concatenate(([0],li))
        plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])
    plot(0.3*np.arange(1000)**(-4/3.), 'k--')
    for t, c in zip([5,13,148], 'bgr'): 
        axvline(t, color=c)
    xscale('log')
    yscale('log')
    ylim(1e-4,1e-1)
    xlim(1,1e3)
    xlabel(r'$t/\tau_B$')
    ylabel(r'strand rupture per particle')
figure(figsize=(12,12))
for i,L in enumerate([3,4,8,16]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x in zip(hists, t0s, xs)
        ]

    #153A
    rates[4] = np.vstack((rates[0], rates[-2]*[1,1]+[300-64,0], rates[-1]*[3,1]+[467-64,0]))
    del rates[-2:]

    #163A
    rates[0] = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    del rates[-1]

    #172A
    rates[2] = np.vstack((rates[2], rates[-1][:460]*[3,1]+[419-33,0]))
    del rates[-1]

    #162B
    rates[1] = np.vstack((rates[1], rates[-1]*[1,1]+[35*6-159,0]))
    del rates[-1]

    #150A
    rates[-2] = np.vstack((rates[-2], rates[-1]+[8,0]))
    del rates[-1]
    
    subplot(2,2,i+1)
    for r in rates:
        #10 points per decade
        li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
        li = np.concatenate(([0],li))
        plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])
    plot(0.3*np.arange(1000)**(-4/3.), 'k--')
    for t, c in zip([5,13,148], 'bgr'): 
        axvline(t, color=c)
    xscale('log')
    yscale('log')
    ylim(1e-5,1e-1)
    xlim(1,1e3)
    xlabel(r'$t/\tau_B$')
    ylabel(r'strand rupture per particle')
t0 = 64
his = hists[4]
plot(np.arange(len(his)-t0), his[t0:,4:].sum(-1))
his = hists[10]
plot(np.arange(len(his))+300-64, his[:,4:].sum(-1))
his = hists[11]
plot(np.arange(len(his))+50*6-64, his[:,4:].sum(-1)/12)
t0 = 64
his = hists[4]
plot(np.arange(len(his)-t0), his[t0:,4:].sum(-1))
his = hists[10]
plot(np.arange(len(his))+300-64, his[:,4:].sum(-1))
his = hists[11]
plot(np.arange(len(his))+110*6-64, his[:,4:].sum(-1)/12)
t0 = 64
his = hists[4]
plot(np.arange(len(his)-t0), his[t0:,4:].sum(-1))
his = hists[10]
plot(np.arange(len(his))+300-64, his[:,4:].sum(-1))
his = hists[11]
plot(np.arange(len(his))+300+167-64, his[:,4:].sum(-1)/12)
t0s = [68,151,33,64,200,104] + [0]*len(hists)
rates = [
    np.column_stack((
        np.arange(len(his)-t0),
        his[t0:,4:].sum(-1)*1./x.get_nb()[t0:-1]
        )) 
    for his, t0, x in zip(hists, t0s, xs)
    ]
print len(rates)
#153A
rates[4] = np.vstack((rates[0], rates[-2]*[1,1]+[300-64,0], rates[-1]*[3,1]+[467-64,0]))
del rates[-2:]
print len(rates)
#163A
rates[0] = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
del rates[-1]
print len(rates)
#172A
rates[2] = np.vstack((rates[2], rates[-1][:460]*[3,1]+[419-33,0]))
del rates[-1]
print len(rates)
#162B
rates[1] = np.vstack((rates[1], rates[-1]*[1,1]+[35*6-159,0]))
del rates[-1]
print len(rates)
#150A
rates[-2] = np.vstack((rates[-2], rates[-1]+[8,0]))
del rates[-1]
print len(rates)
for r in rates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
#yscale('log')
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle.pdf')
for r,nam in zip(rates, names):
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])
    np.savetxt(
        os.path.join('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors', name+'_log.broken'),
        np.column_stack((
            r[li[:-1],0],
            [r[i:j,1].mean() for i,j in zip(li, li[1:])]
            ))
        )
plot(0.3*np.arange(1000)**(-4/3.), 'k--')
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
yscale('log')
ylim(3e-4,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
pa = '/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data'
ns = [
    '3_DenseGel/163A/1_percolation/1340_percolation',
    '2_MidGel/162B/1_percolation/1745',
    '1_DiluteGel/172A/1_percolation/172A_1206_percolation',
    #'1_DiluteGel/155C/1_percolation/155C_percolation_1645',
    '2_MidGel/170A/1_percolation/170A_1451_percolation',
    '2_MidGel/153A/1_percolation/153A_percolation_1540',
    '3_DenseGel/150A/1_percolation/150A_percolation_1227',
    '3_DenseGel/150A/1_percolation/150A_percolation_1239',
    '3_DenseGel/150A/2_ageing/150A_Ageing_1303',
    '2_MidGel/162B/2_ageing/1815_ageing',
    '1_DiluteGel/172A/2_ageing/172A_1318_ageing',
    '3_DenseGel/163A/2_ageing/1415_ageing',
    '2_MidGel/153A/2_percolation/153A_percolation_1610',
    '2_MidGel/153A/3_ageing/153A_ageing_1630',
    ]
xs = [
    xp.Experiment(os.path.join(pa, n+'.traj'))
    for n in ns
    ]
hists =[]
for x in xs:
    lengths = [
        np.loadtxt(name, int) 
        for t, name in x.enum('_broken', 'length') 
        if t>0]
    maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
    hists.append(np.array([
        np.histogram(l, np.arange(-1, maxlength+1))[0]
        for l in lengths]))
len(hists[5]), len(hists[5])+len(hists[6])
t0 = 104
his = hists[5]
plot(np.arange(len(his)-t0), his[t0:,4:].sum(-1))
his = hists[6]
plot(np.arange(len(his))+112-104, his[:,4:].sum(-1))
his = hists[7]
plot(np.arange(len(his))+328-104, his[:,4:].sum(-1)/12)
t0s = [68,151,33,64,200,104] + [0]*len(hists)
rates = [
    np.column_stack((
        np.arange(len(his)-t0),
        his[t0:,4:].sum(-1)*1./x.get_nb()[t0:-1]
        )) 
    for his, t0, x in zip(hists, t0s, xs)
    ]
print len(rates)
#153A
rates[4] = np.vstack((rates[0], rates[-2]*[1,1]+[300-64,0], rates[-1]*[3,1]+[467-64,0]))
del rates[-2:]
print len(rates)
#163A
rates[0] = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
del rates[-1]
print len(rates)
#172A
rates[2] = np.vstack((rates[2], rates[-1][:460]*[3,1]+[419-33,0]))
del rates[-1]
print len(rates)
#162B
rates[1] = np.vstack((rates[1], rates[-1]*[1,1]+[35*6-159,0]))
del rates[-1]
print len(rates)
#150A
rates[5] = np.vstack((rates[5], rates[-2]+[8,0], rates[-2]*[1,1]+[328-104,0]))
del rates[-2:]
print len(rates)
for r in rates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
#yscale('log')
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle.pdf')
for r,nam in zip(rates, names):
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])
    np.savetxt(
        os.path.join('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors', name+'_log.broken'),
        np.column_stack((
            r[li[:-1],0],
            [r[i:j,1].mean() for i,j in zip(li, li[1:])]
            ))
        )
plot(0.3*np.arange(1000)**(-4/3.), 'k--')
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
yscale('log')
ylim(3e-4,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
pa = '/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data'
ns = [
    '3_DenseGel/163A/1_percolation/1340_percolation',
    '2_MidGel/162B/1_percolation/1745',
    '1_DiluteGel/172A/1_percolation/172A_1206_percolation',
    #'1_DiluteGel/155C/1_percolation/155C_percolation_1645',
    '2_MidGel/170A/1_percolation/170A_1451_percolation',
    '2_MidGel/153A/1_percolation/153A_percolation_1540',
    '3_DenseGel/150A/1_percolation/150A_percolation_1227',
    '3_DenseGel/150A/1_percolation/150A_percolation_1239',
    #'3_DenseGel/150A/2_ageing/150A_Ageing_1303',
    '2_MidGel/162B/2_ageing/1815_ageing',
    '1_DiluteGel/172A/2_ageing/172A_1318_ageing',
    '3_DenseGel/163A/2_ageing/1415_ageing',
    '2_MidGel/153A/2_percolation/153A_percolation_1610',
    '2_MidGel/153A/3_ageing/153A_ageing_1630',
    ]
xs = [
    xp.Experiment(os.path.join(pa, n+'.traj'))
    for n in ns
    ]
hists =[]
for x in xs:
    lengths = [
        np.loadtxt(name, int) 
        for t, name in x.enum('_broken', 'length') 
        if t>0]
    maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
    hists.append(np.array([
        np.histogram(l, np.arange(-1, maxlength+1))[0]
        for l in lengths]))
t0s = [68,151,33,64,200,104] + [0]*len(hists)
rates = [
    np.column_stack((
        np.arange(len(his)-t0),
        his[t0:,4:].sum(-1)*1./x.get_nb()[t0:-1]
        )) 
    for his, t0, x in zip(hists, t0s, xs)
    ]
print len(rates)
#153A
rates[4] = np.vstack((rates[0], rates[-2]*[1,1]+[300-64,0], rates[-1]*[3,1]+[467-64,0]))
del rates[-2:]
print len(rates)
#163A
rates[0] = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
del rates[-1]
print len(rates)
#172A
rates[2] = np.vstack((rates[2], rates[-1][:460]*[3,1]+[419-33,0]))
del rates[-1]
print len(rates)
#162B
rates[1] = np.vstack((rates[1], rates[-1]*[1,1]+[35*6-159,0]))
del rates[-1]
print len(rates)
#150A
rates[-2] = np.vstack((rates[-2], rates[-1]+[8,0]))
del rates[-1]
#rates[5] = np.vstack((rates[5], rates[-2]+[8,0], rates[-2]*[1,1]+[328-104,0]))
#del rates[-2:]
print len(rates)
for r in rates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
#yscale('log')
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle.pdf')
for r,nam in zip(rates, names):
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])
    np.savetxt(
        os.path.join('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors', name+'_log.broken'),
        np.column_stack((
            r[li[:-1],0],
            [r[i:j,1].mean() for i,j in zip(li, li[1:])]
            ))
        )
plot(0.3*np.arange(1000)**(-4/3.), 'k--')
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
yscale('log')
ylim(3e-4,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
for r,name in zip(rates, names):
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])
    np.savetxt(
        os.path.join('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors', name+'_log.broken'),
        np.column_stack((
            r[li[:-1],0],
            [r[i:j,1].mean() for i,j in zip(li, li[1:])]
            ))
        )
plot(0.3*np.arange(1000)**(-4/3.), 'k--')
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
yscale('log')
ylim(3e-4,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
get_ipython().system(u'rm /media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors/*_log.broken')
for r,name in zip(rates, names):
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])
    np.savetxt(
        os.path.join('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors', name[:4]+'_log.broken'),
        np.column_stack((
            r[li[:-1],0],
            [r[i:j,1].mean() for i,j in zip(li, li[1:])]
            ))
        )
plot(0.3*np.arange(1000)**(-4/3.), 'k--')
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
yscale('log')
ylim(3e-4,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
for r,name in zip(rates, names):
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])
    np.savetxt(
        os.path.join('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors', name[:4]+'_log.broken'),
        np.column_stack((
            r[li[:-1],0],
            [r[i:j,1].mean() for i,j in zip(li, li[1:])]
            ))
        )
plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(0.3*np.arange(1000)**(-1.57/3.-1), 'k:')
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
yscale('log')
ylim(3e-4,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
for r,name in zip(rates, names):
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])
    np.savetxt(
        os.path.join('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors', name[:4]+'_log.broken'),
        np.column_stack((
            r[li[:-1],0],
            [r[i:j,1].mean() for i,j in zip(li, li[1:])]
            ))
        )
plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(0.4*np.arange(1000)**(-1.57/3.-1), 'k:')
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
yscale('log')
ylim(3e-4,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
for r,name in zip(rates, names):
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])
    np.savetxt(
        os.path.join('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors', name[:4]+'_log.broken'),
        np.column_stack((
            r[li[:-1],0],
            [r[i:j,1].mean() for i,j in zip(li, li[1:])]
            ))
        )
plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(0.5*np.arange(1000)**(-1.57/3.-1), 'k:')
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
yscale('log')
ylim(3e-4,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
t0s = [68,151,33,64,200,104] + [0]*len(hists)
rates = [
    np.column_stack((
        np.arange(len(his)-t0),
        his[t0:,4:].sum(-1)*1./x.get_nb()[t0:-1]
        )) 
    for his, t0, x in zip(hists, t0s, xs)
    ]
print len(rates)
#153A
rates[4] = np.vstack((rates[4], rates[-2]*[1,1]+[300-64,0], rates[-1]*[3,1]+[467-64,0]))
del rates[-2:]
print len(rates)
#163A
rates[0] = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
del rates[-1]
print len(rates)
#172A
rates[2] = np.vstack((rates[2], rates[-1][:460]*[3,1]+[419-33,0]))
del rates[-1]
print len(rates)
#162B
rates[1] = np.vstack((rates[1], rates[-1]*[1,1]+[35*6-159,0]))
del rates[-1]
print len(rates)
#150A
rates[-2] = np.vstack((rates[-2], rates[-1]+[8,0]))
del rates[-1]
#rates[5] = np.vstack((rates[5], rates[-2]+[8,0], rates[-2]*[1,1]+[328-104,0]))
#del rates[-2:]
print len(rates)
for r in rates:
    plot(r[:,0], r[:,1])

for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
#yscale('log')
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle.pdf')
for r,name in zip(rates, names):
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])
    np.savetxt(
        os.path.join('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors', name[:4]+'_log.broken'),
        np.column_stack((
            r[li[:-1],0],
            [r[i:j,1].mean() for i,j in zip(li, li[1:])]
            ))
        )
plot(0.3*np.arange(1000)**(-4/3.), 'k--')
plot(0.5*np.arange(1000)**(-1.57/3.-1), 'k:')
for t, c in zip([5,13,148], 'bgr'): 
    axvline(t, color=c)
xscale('log')
yscale('log')
ylim(3e-4,0.025)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
qmaxs = []
names = [
    '163A_1340_percolation',
    '162B_percolation',
    '172A_1206_percolation',
    '170A_1451_percolation',
    '153A_percolation_1540',
    '150A_percolation_1227',
    '150A_percolation_1239',
    '162B_1815_ageing',
    '172A_1318_ageing',
    '163A_1415_ageing',
    '153A_percolation_1610',
    '153A_ageing_1630',
    ]
imax = 44
for name, t0, x in zip(names, t0s, xs):
    Ss = np.load(os.path.join(
        '/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors',
        name+'_structure_factor.npy'
        ))
    q = np.arange(Ss.shape[1])*Ss[0,0]*2*x.radius
    qmaxs.append(
        np.column_stack((
            np.arange(len(Ss)-t0),
            (Ss[t0:,4:imax]*q[4:imax]).sum(1)/Ss[t0:,4:imax].sum(1)
            ))
        )
#153A
qmaxs[4] = np.vstack((qmaxs[4], qmaxs[-2]*[1,1/2.]+[300-64,0], qmaxs[-1]*[3,1]+[467-64,0]))
del qmaxs[-2:]
#172A
qmaxs[0] = np.vstack((qmaxs[0], qmaxs[-1]*[3,1]+[205-68,0]))
del qmaxs[-1]
#172A
qmaxs[2] = np.vstack((qmaxs[2], qmaxs[-1][:460]*[3,1]+[419-33,0]))
del qmaxs[-1]
#162B
qmaxs[1] = np.vstack((qmaxs[1], qmaxs[-1]*[1,1]+[35*6-159,0]))
del qmaxs[-1]
#150A
qmaxs[5] = np.vstack((qmaxs[5], qmaxs[6]+[8,0]))
del qmaxs[6]
qmaxs = []
names = [
    '163A_1340_percolation',
    '162B_percolation',
    '172A_1206_percolation',
    '170A_1451_percolation',
    '153A_percolation_1540',
    '150A_percolation_1227',
    '150A_percolation_1239',
    '162B_1815_ageing',
    '172A_1318_ageing',
    '163A_1415_ageing',
    '153A_percolation_1610',
    '153A_ageing_1630',
    ]
imax = 44
for name, t0, x in zip(names, t0s, xs):
    Ss = np.load(os.path.join(
        '/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/structure_factors',
        name+'_structure_factor.npy'
        ))
    q = np.arange(Ss.shape[1])*Ss[0,0]*2*x.radius
    qmaxs.append(
        np.column_stack((
            np.arange(len(Ss)-t0),
            (Ss[t0:,4:imax]*q[4:imax]).sum(1)/Ss[t0:,4:imax].sum(1)
            ))
        )
#153A
qmaxs[4] = np.vstack((qmaxs[4], qmaxs[-2]*[1,1]+[300-64,0], qmaxs[-1]*[3,1]+[467-64,0]))
del qmaxs[-2:]
#172A
qmaxs[0] = np.vstack((qmaxs[0], qmaxs[-1]*[3,1]+[205-68,0]))
del qmaxs[-1]
#172A
qmaxs[2] = np.vstack((qmaxs[2], qmaxs[-1][:460]*[3,1]+[419-33,0]))
del qmaxs[-1]
#162B
qmaxs[1] = np.vstack((qmaxs[1], qmaxs[-1]*[1,1]+[35*6-159,0]))
del qmaxs[-1]
#150A
qmaxs[5] = np.vstack((qmaxs[5], qmaxs[6]+[8,0]))
del qmaxs[6]
for qmax in qmaxs:
    plot(qmax[:,0], qmax[:,1])
#plot(np.arange(len(qmaxs[-1]))+8, qmaxs[-1], gca().lines[-1].get_color())
xscale('log');
yscale('log');
xlabel(r'$t/\tau_B$')
ylabel(r'$q_\mathrm{max}\sigma$')
figure(figsize=(12,12))
for i,L in enumerate([3,4,8,16]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x in zip(hists, t0s, xs)
        ]

    #153A
    rates[4] = np.vstack((rates[0], rates[-2]*[1,1]+[300-64,0], rates[-1]*[3,1]+[467-64,0]))
    del rates[-2:]

    #163A
    rates[0] = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    del rates[-1]

    #172A
    rates[2] = np.vstack((rates[2], rates[-1][:460]*[3,1]+[419-33,0]))
    del rates[-1]

    #162B
    rates[1] = np.vstack((rates[1], rates[-1]*[1,1]+[35*6-159,0]))
    del rates[-1]

    #150A
    rates[-2] = np.vstack((rates[-2], rates[-1]+[8,0]))
    del rates[-1]
    
    subplot(2,2,i+1)
    for r in rates:
        #10 points per decade
        li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
        li = np.concatenate(([0],li))
        plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])
    plot(0.3*np.arange(1000)**(-4/3.), 'k--')
    for t, c in zip([5,13,148], 'bgr'): 
        axvline(t, color=c)
    xscale('log')
    yscale('log')
    ylim(1e-5,1e-1)
    xlim(1,1e3)
    xlabel(r'$t/\tau_B$')
    ylabel(r'strand rupture per particle')
for i,L in enumerate([3,4,8,16]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])
for i,L in enumerate([3,4,8,16]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])
xscale('log')
yscale('log')
#ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
for i,L in enumerate([1,2,3,4,8,16]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])
xscale('log')
yscale('log')
#ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
for i,L in enumerate([0,1,2,3,4,8,16]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])
xscale('log')
yscale('log')
#ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
for i,L in enumerate([0,1,3,4,8,16]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])
xscale('log')
yscale('log')
#ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
for i,L in enumerate([0,2,3,4,8,16]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])
xscale('log')
yscale('log')
#ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
for i,L in enumerate([0,1,3,4,8,16]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])
xscale('log')
yscale('log')
#ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
for i,L in enumerate([0,1,3,4,8,16]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(i/5.))
xscale('log')
yscale('log')
#ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
for u,L in enumerate([0,1,3,4,8,16]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/5.))
xscale('log')
yscale('log')
#ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
his.shape
for u,L in enumerate([0,1,3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/5.))
xscale('log')
yscale('log')
#ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
for u,L in enumerate([0,1,3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/6.))
xscale('log')
yscale('log')
#ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
plot(hists[0][:, 0])
plot(hists[0][:, 0])
xscale('log')
yscale('log')
plot(hists[0][:, 0])
plot(hists[0][:, 1])
xscale('log')
yscale('log')
plot(hists[0][:, 0])
plot(hists[0][:, 1])
plot(hists[0][:, 2])
xscale('log')
yscale('log')
plot(hists[0][:, 0])
plot(hists[0][:, 1])
plot(hists[0][:, 3])
xscale('log')
yscale('log')
plot(hists[0][:, 0])
plot(hists[0][:, 1])
plot(hists[0][:, 3])
plot(hists[0][:, 4:])
xscale('log')
yscale('log')
plot(hists[0][:, 0])
plot(hists[0][:, 1])
plot(hists[0][:, 3])
plot(hists[0][:, 4:].sum(-1))
xscale('log')
yscale('log')
plot(hists[0][:, 0]*1./x.get_nb()[:-1])
plot(hists[0][:, 1]*1./x.get_nb()[:-1])
plot(hists[0][:, 3]*1./x.get_nb()[:-1])
plot(hists[0][:, 4:].sum(-1)*1./x.get_nb()[:-1])
xscale('log')
yscale('log')
plot(hists[0][:, 0]*1./xs[0].get_nb()[:-1])
plot(hists[0][:, 1]*1./xs[0].get_nb()[:-1])
plot(hists[0][:, 3]*1./xs[0].get_nb()[:-1])
plot(hists[0][:, 4:].sum(-1)*1./xs[0].get_nb()[:-1])
xscale('log')
yscale('log')
for h in [
    hists[0][:, 0],
    hists[0][:, 1],
    hists[0][:, 3],
    hists[0][:, 4:].sum(-1)
    ]:
    plot(h*1./xs[0].get_nb()[:-1])
xscale('log')
yscale('log')
for h in [
    hists[0][:, 0],
    hists[0][:, 1],
    hists[0][:, 3],
    hists[0][:, 4:].sum(-1)
    ]:
    r = h*1./xs[0].get_nb()[:-1]
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(li[:-1], [r[i:j].mean() for i,j in zip(li, li[1:])])
xscale('log')
yscale('log')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/6.))
xscale('log')
yscale('log')
#ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.))
xscale('log')
yscale('log')
#ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.))

plot(0.3*np.arange(1000)**(-1.5), 'k--')
xscale('log')
yscale('log')
#ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.))

plot(0.3*np.arange(1000)**(-1.5), 'k--')
plot(0.3*np.arange(1000)**(-2), 'k--')
xscale('log')
yscale('log')
#ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.))

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-2), 'k--')
xscale('log')
yscale('log')
#ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.))

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-2), 'k:')
xscale('log')
yscale('log')
#ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.))

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 0.1*np.arange(1,1000)**(-2), 'k:')
xscale('log')
yscale('log')
#ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.))

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 0.1*np.arange(1,1000)**(-2.), 'k:')
xscale('log')
yscale('log')
#ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.))

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 0.1*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
#ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.))

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
#ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.))

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.))

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-4/3.), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.))
    gca().annotate(
        'L>%d'%(i-3), 
        xy=gca().lines[0].get_data()[-1], 
        xytext=(3, 1.5),
        )

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-4/3.), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.))
    gca().annotate(
        'L>%d'%(i-3), 
        xy=gca().lines[0].get_data()[-1], 
        xytext=(3, 1.5),
        )

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-4/3.), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
#savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.))
    gca().annotate(
        'L>%d'%(i-3), 
        xy=gca().lines[0].get_data()[-1], xycoords='data',
        xytext=(0, 0),
        )

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-4/3.), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
#savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.))
    gca().annotate(
        'L>%d'%(i-3), 
        xy=gca().lines[0].get_data()[-1], xycoords='data',
        xytext=(0, 0),
        )

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-4/3.), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
get_ipython().magic(u'debug')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.))
    gca().annotate(
        'L>%d'%(i-3), 
        xy=gca().lines[0].get_data()[:,-1], xycoords='data',
        xytext=(0, 0),
        )

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-4/3.), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
get_ipython().magic(u'debug')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.))
    print l.shape
    gca().annotate(
        'L>%d'%(i-3), 
        xy=l.get_data()[:,-1], xycoords='data',
        xytext=(0, 0),
        )

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-4/3.), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.))
    print l.shape
    gca().annotate(
        'L>%d'%(i-3), 
        xy=(l.get_xdata()[-1], l.get_ydata()[-1]), xycoords='data',
        xytext=(0, 0),
        )

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-4/3.), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.))
    gca().annotate(
        'L>%d'%(i-3), 
        xy=(l.get_xdata()[-1], l.get_ydata()[-1]), xycoords='data',
        xytext=(0, 0),
        )

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-4/3.), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.))[0]
    gca().annotate(
        'L>%d'%(i-3), 
        xy=(l.get_xdata()[-1], l.get_ydata()[-1]), xycoords='data',
        xytext=(0, 0),
        )

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-4/3.), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.))[0]
    print l.get_data()
    print l.get_xdata()
    print l.get_ydata()
    gca().annotate(
        'L>%d'%(i-3), 
        xy=(l.get_xdata()[-1], l.get_ydata()[-1]), xycoords='data',
        xytext=(0, 0),
        )

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-4/3.), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
for u,L in enumerate([3,4,8,16]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.))[0]
    print l.get_data()
    print l.get_xdata()
    print l.get_ydata()
    gca().annotate(
        'L>%d'%(i-3), 
        xy=(l.get_xdata()[-1], l.get_ydata()[-1]), xycoords='data',
        xytext=(0, 0),
        )

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-4/3.), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.))[0]
    print l.get_data()
    print l.get_xdata()
    print l.get_ydata()
    #gca().annotate(
     #   'L>%d'%(i-3), 
      #  xy=(l.get_xdata()[-1], l.get_ydata()[-1]), xycoords='data',
       # xytext=(0, 0),
        #)

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-4/3.), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.))[0]
    end = np.where(l.get_ydata()>0)[0][-1]
    gca().annotate(
        'L>%d'%(i-3), 
        xy=(l.get_xdata()[end], l.get_ydata()[end]), xycoords='data',
        xytext=(0, 0),
        )

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-4/3.), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.))[0]
    end = np.where(l.get_ydata()>0)[0][-1]
    gca().annotate(
        'L>%d'%(i-3), 
        xy=(l.get_xdata()[end], l.get_ydata()[end]), xycoords='data',
        xytext=(0, 0),
        )

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-4/3.), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
#savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.), label='L>%d'%(L-3))

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-4/3.), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
#savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.), label='L>%d'%(L-3))

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-4/3.), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
legend(log='upper right')
#savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.), label='L>%d'%(L-3))

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-4/3.), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
legend(loc='upper right')
#savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.), label='L>%d'%(L-3))

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-4/3.), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
legend(loc='upper right', ncols=5)
#savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.), label='L>%d'%(L-3))

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-4/3.), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
legend(loc='upper right', ncol=5)
#savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.))[0]
    end = np.where(l.get_ydata()>0)[0][-1]
    annotate('L>%d'%(L-3), xy=(l.get_xdata()[end], l.get_ydata()[end]))

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-4/3.), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
#savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.))[0]
    end = np.where(l.get_ydata()>0)[0][-1]
    text(l.get_xdata()[end], l.get_ydata()[end], 'L>%d'%(L-3))

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-4/3.), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
#savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.))[0]
    end = np.where(l.get_ydata()>0)[0][-1]
    text(l.get_xdata()[end], l.get_ydata()[end], 'L>%d'%(L-3), color=cm.autumn(u/4.))

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-4/3.), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
#savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.))[0]
    end = np.where(l.get_ydata()>0)[0][-1]
    text(1.1*l.get_xdata()[end], l.get_ydata()[end], 'L>%d'%(L-3), color=cm.autumn(u/4.))

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-4/3.), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
#savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.))[0]
    end = np.where(l.get_ydata()>0)[0][-1]
    text(1.1*l.get_xdata()[end], l.get_ydata()[end], 'L>%d'%(L-3), color=cm.autumn(u/4.))

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
#savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.))[0]
    end = np.where(l.get_ydata()>0)[0][-1]
    text(1.1*l.get_xdata()[end], l.get_ydata()[end], 'L>%d'%(L-3), color=cm.autumn(u/4.), verticalalignement='top')

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
#savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.))[0]
    end = np.where(l.get_ydata()>0)[0][-1]
    text(1.1*l.get_xdata()[end], l.get_ydata()[end], 'L>%d'%(L-3), color=cm.autumn(u/4.), verticalalignment='top')

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
#savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.))[0]
    end = np.where(l.get_ydata()>0)[0][-1]
    txt = text(1.1*l.get_xdata()[end], l.get_ydata()[end], 'L>%d'%(L-3), color=cm.autumn(u/4.))

txt.set_color('k')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
#savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.))[0]
    end = np.where(l.get_ydata()>0)[0][-1]
    txt = text(1.1*l.get_xdata()[end], l.get_ydata()[end], 'L>%d'%(L-3), color=cm.autumn(u/4.))

txt.set_color('k')
txt.set_horizontalalignment('center')
txt.set_verticalalalignment('top')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
#savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.))[0]
    end = np.where(l.get_ydata()>0)[0][-1]
    txt = text(1.1*l.get_xdata()[end], l.get_ydata()[end], 'L>%d'%(L-3), color=cm.autumn(u/4.))

txt.set_color('k')
txt.set_horizontalalignment('center')
txt.set_verticalalignment('top')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
#savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
for u,L in enumerate([3,4,8,16,32]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.))[0]
    end = np.where(l.get_ydata()>0)[0][-1]
    txt = text(1.1*l.get_xdata()[end], l.get_ydata()[end], 'L>%d'%(L-3), color=cm.autumn(u/4.))

txt.set_color('k')
txt.set_horizontalalignment('center')
txt.set_verticalalignment('top')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
for u,L in enumerate([3,4,5,8,13,23,33]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/4.))[0]
    end = np.where(l.get_ydata()>0)[0][-1]
    txt = text(1.1*l.get_xdata()[end], l.get_ydata()[end], 'L>%d'%(L-3), color=cm.autumn(u/4.))

txt.set_color('k')
txt.set_horizontalalignment('center')
txt.set_verticalalignment('top')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
ylim(1e-5,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
for u,L in enumerate([3,4,5,8,13,23,33]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/7.))[0]
    end = np.where(l.get_ydata()>0)[0][-1]
    txt = text(1.1*l.get_xdata()[end], l.get_ydata()[end], 'L>%d'%(L-3), color=cm.autumn(u/7.))

txt.set_color('k')
txt.set_horizontalalignment('center')
txt.set_verticalalignment('top')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
ylim(1e-6,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
for u,L in enumerate([3,4,5,8,13,23,33]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/6.))[0]
    end = np.where(l.get_ydata()>0)[0][-1]
    txt = text(1.1*l.get_xdata()[end], l.get_ydata()[end], 'L>%d'%(L-3), color=cm.autumn(u/6.))

txt.set_color('k')
txt.set_horizontalalignment('center')
txt.set_verticalalignment('top')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
ylim(1e-6,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
for u,L in enumerate([3,4,5,8,13,23,33]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/6.))[0]
    end = np.where(l.get_ydata()>0)[0][-1]
    txt = text(1.1*l.get_xdata()[end], l.get_ydata()[end], 'L>%d'%(L-3), color=cm.autumn(u/6.))

txt.set_color('k')
txt.set_horizontalalignment('left')
txt.set_verticalalignment('top')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
ylim(1e-6,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
for u,L in enumerate([3,4,5,8,13,23,33]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/6.))[0]
    end = np.where(l.get_ydata()>0)[0][-1]
    txt = text(1.1*l.get_xdata()[end], l.get_ydata()[end], 'L>%d'%(L-3), color=cm.autumn(u/6.))

txt.set_color('k')
txt.set_horizontalalignment('right')
txt.set_verticalalignment('top')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
ylim(1e-6,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
for u,L in enumerate([3,4,5,8,13,23,33]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/6.))[0]
    end = np.where(l.get_ydata()>0)[0][-1]
    txt = text(1.1*l.get_xdata()[end], l.get_ydata()[end], 'L>%d'%(L-3), color=cm.autumn(u/6.))

txt.set_color('k')
txt.set_horizontalalignment('right')
txt.set_verticalalignment('top')
for a in np.arange(1.5,4,0.5):
    plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-a), 'k--')
xscale('log')
yscale('log')
ylim(1e-6,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
for u,L in enumerate([3,4,5,8,13,23,33]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/6.))[0]
    end = np.where(l.get_ydata()>0)[0][-1]
    txt = text(1.1*l.get_xdata()[end], l.get_ydata()[end], 'L>%d'%(L-3), color=cm.autumn(u/6.))

txt.set_color('k')
txt.set_horizontalalignment('right')
txt.set_verticalalignment('top')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-3.), 'k:')
xscale('log')
yscale('log')
ylim(1e-6,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
for u,L in enumerate([3,4,5,8,13,23,33]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/6.))[0]
    end = np.where(l.get_ydata()>0)[0][-1]
    txt = text(1.1*l.get_xdata()[end], l.get_ydata()[end], 'L>%d'%(L-3), color=cm.autumn(u/6.))

txt.set_color('k')
txt.set_horizontalalignment('right')
txt.set_verticalalignment('top')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-5.), 'k:')
xscale('log')
yscale('log')
ylim(1e-6,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
for u,L in enumerate([3,4,5,8,13,23,33]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/6.))[0]
    end = np.where(l.get_ydata()>0)[0][-1]
    txt = text(1.1*l.get_xdata()[end], l.get_ydata()[end], 'L>%d'%(L-3), color=cm.autumn(u/6.))

txt.set_color('k')
txt.set_horizontalalignment('right')
txt.set_verticalalignment('top')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 3*np.arange(1,1000)**(-5.), 'k:')
xscale('log')
yscale('log')
ylim(1e-6,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
for u,L in enumerate([3,4,5,8,13,23,33]):
    rates = [
        np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )) 
        for his, t0, x, name in zip(hists, t0s, xs, names)
        if name[:4] == '163A'
        ]
    #concatenate
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/6.))[0]
    end = np.where(l.get_ydata()>0)[0][-1]
    txt = text(1.1*l.get_xdata()[end], l.get_ydata()[end], 'L>%d'%(L-3), color=cm.autumn(u/6.))

txt.set_color('k')
txt.set_horizontalalignment('right')
txt.set_verticalalignment('top')
plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 10*np.arange(1,1000)**(-5.), 'k:')
xscale('log')
yscale('log')
ylim(1e-6,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/breaking_per_particle_Lthr.pdf')
import networkx as nx
def shortest_length(g, trs):
    """From a graph g and couples of nodes trs, return the on-graph distance between each couple. 
    Returns -1 if either member of the couple do not belong to the graph.
    
    trs is a (N,2) array of integers."""
    lengths = np.zeros(len(trs))
    for i, (a,b) in enumerate(trs):
        if not g.has_node(a) or not g.has_node(b): continue
        try:
            lengths[i] = nx.shortest_path_length(g, a,b)
        except nx.NetworkXNoPath:
            lengths[i] = -1
    return lengths
def broken_bonds_lenghts(bonds0, bonds1, p2tr0, p2tr1):
    """On graph distance at t1 between particles involved in a bond at t0 and no more at t1.
    
    bonds0, bonds1 are respectively the bonds at t0 and 1 in terms of position
    p2tr0, p2tr1 are respectively the position to trajectory relationship at t0 and t1"""
    #bonds (between trajectories) existing at t but no more at t+dt
    # = broken bonds + lost trajectories
    trbonds = set([(a,b) for a,b in np.sort(p2tr0[bonds0], 1)]) - set([(a,b) for a,b in np.sort(p2tr1[bonds1], axis=1)])

    #graph of the bonds between trajectories at t+dt
    g = nx.Graph()
    g.add_nodes_from(p2tr1)
    g.add_edges_from(p2tr1[bonds1])

    return shortest_length(g, trbonds)
x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
m = np.ones(3)*1024
M = np.zeros(3)
for t, name in x.enum():
    pos = np.loadtxt(name, skiprows=2)
    m = np.minimum(m, pos.min(0))
    M = np.maximum(M, pos.max(0))
import shutil
x.tx.trajfile
x.trajfile
LL = 256
newpath = os.path.join(x.path, '%d'%L)
os.mkdir(newpath)
shutil.copy(os.path.join(x.path, x.trajfile), newpath)
for t, name in x.enum(ext='p2tr'):
    shutil.copy(name, newpath)

c = 0.5*(m+M)
for t, name in x.enum(ext='bonds'):
    pos = np.loadtxt(x.get_format_string()%t, skiprows=2)
    inside = np.min(np.abs(pos-c)<LL, -1)
    try:
        bonds = np.atleast_2d(np.loadtxt(name, dtype=int))
    except UserWarning:
        bonds = np.zeros([0,2], int)
    #select bonds that are inside
    bonds = bonds[inside[bonds].min(-1)]
    np.savetxt(
        os.path.join(newpath, os.path.split(name)[-1]),
        bonds,
        fmt='%d'
        )
LL = 256
newpath = os.path.join(x.path, '%d'%LL)
os.mkdir(newpath)
os.symlink(os.path.join(x.path, x.trajfile), os.path.join(newpath, x.trajfile))
for t, name in x.enum(ext='p2tr'):
    os.symlink(name, os.path.join(newpath, os.path.split(name)[-1]))

c = 0.5*(m+M)
for t, name in x.enum(ext='bonds'):
    pos = np.loadtxt(x.get_format_string()%t, skiprows=2)
    inside = np.min(np.abs(pos-c)<LL/2, -1)
    try:
        bonds = np.atleast_2d(np.loadtxt(name, dtype=int))
    except UserWarning:
        bonds = np.zeros([0,2], int)
    #select bonds that are inside
    bonds = bonds[inside[bonds].min(-1)]
    np.savetxt(
        os.path.join(newpath, os.path.split(name)[-1]),
        bonds,
        fmt='%d'
        )
x = xp.Experiment(os.path.join(newpath, x.trajfile))
#x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
maxlength = 0
pro = ProgressBar(x.size)

#compute and save the length of broken bonds between consecutive time steps
for t, name in x.enum(ext='bonds'):
    try:
        bonds1 = np.atleast_2d(np.loadtxt(name, dtype=int))
    except UserWarning:
        bonds1 = np.zeros([0,2], int)
    p2tr1 = np.loadtxt(x.get_format_string(ext='p2tr')%t, dtype=int)
    if t>0:
        lengths = broken_bonds_lenghts(bonds0, bonds1, p2tr0, p2tr1)
        np.savetxt(x.get_format_string('_broken', ext='length')%t, lengths, fmt='%d')
        if len(lengths)>0:
            maxlength = max([maxlength, lengths.max()])
    bonds0 = bonds1
    p2tr0 = p2tr1
    pro.animate(t)

#make the time dependent histogram of broken lengths
histlength = np.zeros(maxlength+1, int)
for t, name in x.enum('_broken_q2', ext='length'):
    if t==0: continue
    try:
        lengths = np.loadtxt(name, int)
    except UserWarning:
        continue
    histlength += np.histogram(
        lengths, 
        bins=np.arange(-1, maxlength+1)
        )[0]

#save the time dependent histogram of broken lengths
np.savetxt(
    os.path.join(x.path, 'broken_length.hist'),
    np.column_stack((np.arange(-1, maxlength), histlength)),
    fmt='%d')
x.get_format_string('_broken', ext='length')
histlength = np.zeros(maxlength+1, int)
for t, name in x.enum('_broken', ext='length'):
    if t==0: continue
    try:
        lengths = np.loadtxt(name, int)
    except UserWarning:
        continue
    histlength += np.histogram(
        lengths, 
        bins=np.arange(-1, maxlength+1)
        )[0]

#save the time dependent histogram of broken lengths
np.savetxt(
    os.path.join(x.path, 'broken_length.hist'),
    np.column_stack((np.arange(-1, maxlength), histlength)),
    fmt='%d')
rates = [
    np.column_stack((
        np.arange(len(his)-t0),
        his[t0:,3:].sum(-1)*1./x.get_nb()[t0:-1]
        )) 
    for his, t0, x, name in zip(hists, t0s, [
        xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/%s.traj'%(name)) 
        for name in ['1340_percolation', '256/1340_percolation']
        ])
    ]
for r in rates:
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], color=cm.autumn(u/6.))[0]

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 10*np.arange(1,1000)**(-5.), 'k:')
xscale('log')
yscale('log')
ylim(1e-6,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
rates = [
    np.column_stack((
        np.arange(len(his)-t0),
        his[t0:,3:].sum(-1)*1./x.get_nb()[t0:-1]
        )) 
    for his, t0, x in zip(hists, t0s, [
        xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/%s.traj'%(name)) 
        for name in ['1340_percolation', '256/1340_percolation']
        ])
    ]
for r in rates:
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])[0]

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 10*np.arange(1,1000)**(-5.), 'k:')
xscale('log')
yscale('log')
ylim(1e-6,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
nb = np.zeros(x.size, int)
for t, name in x.enum(ext='bonds'):
    pos = np.loadtxt(x.get_format_string()%t, skiprows=2)
    inside = np.min(np.abs(pos-c)<LL/2, -1)
    nb[t] = inside.sum()
np.savetxt(
    os.path.join(newpath, x.trajfile[:-4]+'nb'),
    nb, fmt='%d'
    )
rates = [
    np.column_stack((
        np.arange(len(his)-t0),
        his[t0:,3:].sum(-1)*1./x.get_nb()[t0:-1]
        )) 
    for his, t0, x in zip(hists, t0s, [
        xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/%s.traj'%name) 
        for name in ['1340_percolation', '256/1340_percolation']
        ])
    ]
for r in rates:
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])[0]

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 10*np.arange(1,1000)**(-5.), 'k:')
xscale('log')
yscale('log')
ylim(1e-6,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
os.path.splitext('dodo.traj')
x.head
reload(xp)
rates = [
    np.column_stack((
        np.arange(len(his)-t0),
        his[t0:,3:].sum(-1)*1./x.get_nb()[t0:-1]
        )) 
    for his, t0, x in zip(hists, t0s, [
        xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/%s.traj'%name) 
        for name in ['1340_percolation', '256/1340_percolation']
        ])
    ]
for r in rates:
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])[0]

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 10*np.arange(1,1000)**(-5.), 'k:')
xscale('log')
yscale('log')
ylim(1e-6,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
rates = [
    np.column_stack((
        np.arange(len(his)-t0),
        his[t0:,3:].sum(-1)*1./x.get_nb()[t0:-1]
        )) 
    for his, t0, x in zip(hists, [0,0], [
        xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/%s.traj'%name) 
        for name in ['1340_percolation', '256/1340_percolation']
        ])
    ]
for r in rates:
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])[0]

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 10*np.arange(1,1000)**(-5.), 'k:')
xscale('log')
yscale('log')
ylim(1e-6,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
for name in ['1340_percolation', '256/1340_percolation']:
    x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/%s.traj'%name)
    lengths = [
        np.loadtxt(name, int) 
        for t, name in x.enum('_broken', 'length') 
        if t>0]
    maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
    hist = np.array([
        np.histogram(l, np.arange(-1, maxlength+1))[0]
        for l in lengths])
    r = his[:,3:].sum(-1)*1./x.get_nb()[:-1]
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(li[:-1], [r[i:j].mean() for i,j in zip(li, li[1:])])[0]

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 10*np.arange(1,1000)**(-5.), 'k:')
xscale('log')
yscale('log')
ylim(1e-6,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
for name in ['1340_percolation', '256/1340_percolation']:
    x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/%s.traj'%name)
    lengths = [
        np.loadtxt(name, int) 
        for t, name in x.enum('_broken', 'length') 
        if t>0]
    maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
    his = np.array([
        np.histogram(l, np.arange(-1, maxlength+1))[0]
        for l in lengths])
    r = his[:,3:].sum(-1)*1./x.get_nb()[:-1]
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(li[:-1], [r[i:j].mean() for i,j in zip(li, li[1:])])[0]

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 10*np.arange(1,1000)**(-5.), 'k:')
xscale('log')
yscale('log')
ylim(1e-6,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
for name, t0 in zip(['1340_percolation', '256/1340_percolation'], [68]*2):
    x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/%s.traj'%name)
    lengths = [
        np.loadtxt(name, int) 
        for t, name in x.enum('_broken', 'length') 
        if t>0]
    maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
    his = np.array([
        np.histogram(l, np.arange(-1, maxlength+1))[0]
        for l in lengths])
    r = his[t0:,3:].sum(-1)*1./x.get_nb()[t0:-1]
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(li[:-1], [r[i:j].mean() for i,j in zip(li, li[1:])])[0]

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 10*np.arange(1,1000)**(-5.), 'k:')
xscale('log')
yscale('log')
ylim(1e-6,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
L = 5
for name, t0 in zip(['1340_percolation', '256/1340_percolation'], [68]*2):
    x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/%s.traj'%name)
    lengths = [
        np.loadtxt(name, int) 
        for t, name in x.enum('_broken', 'length') 
        if t>0]
    maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
    his = np.array([
        np.histogram(l, np.arange(-1, maxlength+1))[0]
        for l in lengths])
    r = his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(li[:-1], [r[i:j].mean() for i,j in zip(li, li[1:])])[0]

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 10*np.arange(1,1000)**(-5.), 'k:')
xscale('log')
yscale('log')
ylim(1e-6,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
m = np.ones(3)*1024
M = np.zeros(3)
for t, name in x.enum():
    pos = np.loadtxt(name, skiprows=2)
    m = np.minimum(m, pos.min(0))
    M = np.maximum(M, pos.max(0))
LL = 128
newpath = os.path.join(x.path, '%d'%LL)
os.mkdir(newpath)
#cheap copy of file.traj
os.symlink(os.path.join(x.path, x.trajfile), os.path.join(newpath, x.trajfile))
#cheap copy of file_t???.p2tr
for t, name in x.enum(ext='p2tr'):
    os.symlink(name, os.path.join(newpath, os.path.split(name)[-1]))

c = 0.5*(m+M)
nb = np.zeros(x.size, int)
#Remove bonds when on of the members are outside
for t, name in x.enum(ext='bonds'):
    pos = np.loadtxt(x.get_format_string()%t, skiprows=2)
    inside = np.min(np.abs(pos-c)<LL/2, -1)
    nb[t] = inside.sum()
    try:
        bonds = np.atleast_2d(np.loadtxt(name, dtype=int))
    except UserWarning:
        bonds = np.zeros([0,2], int)
    #select bonds that are inside
    bonds = bonds[inside[bonds].min(-1)]
    np.savetxt(
        os.path.join(newpath, os.path.split(name)[-1]),
        bonds,
        fmt='%d'
        )
np.savetxt(
    os.path.join(newpath, x.trajfile[:-4]+'nb'),
    nb, fmt='%d'
    )
x = xp.Experiment(os.path.join(newpath, x.trajfile))
#x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
maxlength = 0
pro = ProgressBar(x.size)

#compute and save the length of broken bonds between consecutive time steps
for t, name in x.enum(ext='bonds'):
    try:
        bonds1 = np.atleast_2d(np.loadtxt(name, dtype=int))
    except UserWarning:
        bonds1 = np.zeros([0,2], int)
    p2tr1 = np.loadtxt(x.get_format_string(ext='p2tr')%t, dtype=int)
    if t>0:
        lengths = broken_bonds_lenghts(bonds0, bonds1, p2tr0, p2tr1)
        np.savetxt(x.get_format_string('_broken', ext='length')%t, lengths, fmt='%d')
        if len(lengths)>0:
            maxlength = max([maxlength, lengths.max()])
    bonds0 = bonds1
    p2tr0 = p2tr1
    pro.animate(t)

#make the time dependent histogram of broken lengths
histlength = np.zeros(maxlength+1, int)
for t, name in x.enum('_broken', ext='length'):
    if t==0: continue
    try:
        lengths = np.loadtxt(name, int)
    except UserWarning:
        continue
    histlength += np.histogram(
        lengths, 
        bins=np.arange(-1, maxlength+1)
        )[0]

#save the time dependent histogram of broken lengths
np.savetxt(
    os.path.join(x.path, 'broken_length.hist'),
    np.column_stack((np.arange(-1, maxlength), histlength)),
    fmt='%d')
L = 5
for name, t0 in zip(['1340_percolation', '256/1340_percolation', '128/1340_percolation'], [68]*2):
    x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/%s.traj'%name)
    lengths = [
        np.loadtxt(name, int) 
        for t, name in x.enum('_broken', 'length') 
        if t>0]
    maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
    his = np.array([
        np.histogram(l, np.arange(-1, maxlength+1))[0]
        for l in lengths])
    r = his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(li[:-1], [r[i:j].mean() for i,j in zip(li, li[1:])])[0]

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 10*np.arange(1,1000)**(-5.), 'k:')
xscale('log')
yscale('log')
ylim(1e-6,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
L = 5
for name, t0 in zip(['1340_percolation', '256/1340_percolation', '128/1340_percolation'], [68]*3):
    x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/%s.traj'%name)
    lengths = [
        np.loadtxt(name, int) 
        for t, name in x.enum('_broken', 'length') 
        if t>0]
    maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
    his = np.array([
        np.histogram(l, np.arange(-1, maxlength+1))[0]
        for l in lengths])
    r = his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(li[:-1], [r[i:j].mean() for i,j in zip(li, li[1:])])[0]

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 10*np.arange(1,1000)**(-5.), 'k:')
xscale('log')
yscale('log')
ylim(1e-6,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
LL = 64
newpath = os.path.join(x.path, '%d'%LL)
os.mkdir(newpath)
#cheap copy of file.traj
os.symlink(os.path.join(x.path, x.trajfile), os.path.join(newpath, x.trajfile))
#cheap copy of file_t???.p2tr
for t, name in x.enum(ext='p2tr'):
    os.symlink(name, os.path.join(newpath, os.path.split(name)[-1]))

c = 0.5*(m+M)
nb = np.zeros(x.size, int)
#Remove bonds when on of the members are outside
for t, name in x.enum(ext='bonds'):
    pos = np.loadtxt(x.get_format_string()%t, skiprows=2)
    inside = np.min(np.abs(pos-c)<LL/2, -1)
    nb[t] = inside.sum()
    try:
        bonds = np.atleast_2d(np.loadtxt(name, dtype=int))
    except UserWarning:
        bonds = np.zeros([0,2], int)
    #select bonds that are inside
    bonds = bonds[inside[bonds].min(-1)]
    np.savetxt(
        os.path.join(newpath, os.path.split(name)[-1]),
        bonds,
        fmt='%d'
        )
np.savetxt(
    os.path.join(newpath, x.trajfile[:-4]+'nb'),
    nb, fmt='%d'
    )
#x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
maxlength = 0
pro = ProgressBar(x.size)

#compute and save the length of broken bonds between consecutive time steps
for t, name in x.enum(ext='bonds'):
    try:
        bonds1 = np.atleast_2d(np.loadtxt(name, dtype=int))
    except UserWarning:
        bonds1 = np.zeros([0,2], int)
    p2tr1 = np.loadtxt(x.get_format_string(ext='p2tr')%t, dtype=int)
    if t>0:
        lengths = broken_bonds_lenghts(bonds0, bonds1, p2tr0, p2tr1)
        np.savetxt(x.get_format_string('_broken', ext='length')%t, lengths, fmt='%d')
        if len(lengths)>0:
            maxlength = max([maxlength, lengths.max()])
    bonds0 = bonds1
    p2tr0 = p2tr1
    pro.animate(t)

#make the time dependent histogram of broken lengths
histlength = np.zeros(maxlength+1, int)
for t, name in x.enum('_broken', ext='length'):
    if t==0: continue
    try:
        lengths = np.loadtxt(name, int)
    except UserWarning:
        continue
    histlength += np.histogram(
        lengths, 
        bins=np.arange(-1, maxlength+1)
        )[0]

#save the time dependent histogram of broken lengths
np.savetxt(
    os.path.join(x.path, 'broken_length.hist'),
    np.column_stack((np.arange(-1, maxlength), histlength)),
    fmt='%d')
L = 5
for name, t0 in zip(['1340_percolation', '256/1340_percolation', '128/1340_percolation', '64/1340_percolation'], [68]*3):
    x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/%s.traj'%name)
    lengths = [
        np.loadtxt(name, int) 
        for t, name in x.enum('_broken', 'length') 
        if t>0]
    maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
    his = np.array([
        np.histogram(l, np.arange(-1, maxlength+1))[0]
        for l in lengths])
    r = his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(li[:-1], [r[i:j].mean() for i,j in zip(li, li[1:])])[0]

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 10*np.arange(1,1000)**(-5.), 'k:')
xscale('log')
yscale('log')
ylim(1e-6,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
L = 5
for name, t0 in zip(['1340_percolation', '256/1340_percolation', '128/1340_percolation', '64/1340_percolation'], [68]*4):
    x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/%s.traj'%name)
    lengths = [
        np.loadtxt(name, int) 
        for t, name in x.enum('_broken', 'length') 
        if t>0]
    maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
    his = np.array([
        np.histogram(l, np.arange(-1, maxlength+1))[0]
        for l in lengths])
    r = his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(li[:-1], [r[i:j].mean() for i,j in zip(li, li[1:])])[0]

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 10*np.arange(1,1000)**(-5.), 'k:')
xscale('log')
yscale('log')
ylim(1e-6,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
x = xp.Experiment(os.path.join(newpath, x.trajfile))
#x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
maxlength = 0
pro = ProgressBar(x.size)

#compute and save the length of broken bonds between consecutive time steps
for t, name in x.enum(ext='bonds'):
    try:
        bonds1 = np.atleast_2d(np.loadtxt(name, dtype=int))
    except UserWarning:
        bonds1 = np.zeros([0,2], int)
    p2tr1 = np.loadtxt(x.get_format_string(ext='p2tr')%t, dtype=int)
    if t>0:
        lengths = broken_bonds_lenghts(bonds0, bonds1, p2tr0, p2tr1)
        np.savetxt(x.get_format_string('_broken', ext='length')%t, lengths, fmt='%d')
        if len(lengths)>0:
            maxlength = max([maxlength, lengths.max()])
    bonds0 = bonds1
    p2tr0 = p2tr1
    pro.animate(t)

#make the time dependent histogram of broken lengths
histlength = np.zeros(maxlength+1, int)
for t, name in x.enum('_broken', ext='length'):
    if t==0: continue
    try:
        lengths = np.loadtxt(name, int)
    except UserWarning:
        continue
    histlength += np.histogram(
        lengths, 
        bins=np.arange(-1, maxlength+1)
        )[0]

#save the time dependent histogram of broken lengths
np.savetxt(
    os.path.join(x.path, 'broken_length.hist'),
    np.column_stack((np.arange(-1, maxlength), histlength)),
    fmt='%d')
L = 5
for name, t0 in zip(['1340_percolation', '256/1340_percolation', '128/1340_percolation', '64/1340_percolation'], [68]*4):
    x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/%s.traj'%name)
    lengths = [
        np.loadtxt(name, int) 
        for t, name in x.enum('_broken', 'length') 
        if t>0]
    maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
    his = np.array([
        np.histogram(l, np.arange(-1, maxlength+1))[0]
        for l in lengths])
    r = his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(li[:-1], [r[i:j].mean() for i,j in zip(li, li[1:])])[0]

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 10*np.arange(1,1000)**(-5.), 'k:')
xscale('log')
yscale('log')
ylim(1e-6,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/2_ageing/1415_ageing.traj')
m = np.ones(3)*1024
M = np.zeros(3)
for t, name in x.enum():
    pos = np.loadtxt(name, skiprows=2)
    m = np.minimum(m, pos.min(0))
    M = np.maximum(M, pos.max(0))
x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/2_ageing/1415_ageing.traj')
LL = 256
newpath = os.path.join(x.path, '%d'%LL)
os.mkdir(newpath)
#cheap copy of file.traj
os.symlink(os.path.join(x.path, x.trajfile), os.path.join(newpath, x.trajfile))
#cheap copy of file_t???.p2tr
for t, name in x.enum(ext='p2tr'):
    os.symlink(name, os.path.join(newpath, os.path.split(name)[-1]))

c = 0.5*(m+M)
nb = np.zeros(x.size, int)
#Remove bonds when on of the members are outside
for t, name in x.enum(ext='bonds'):
    pos = np.loadtxt(x.get_format_string()%t, skiprows=2)
    inside = np.min(np.abs(pos-c)<LL/2, -1)
    nb[t] = inside.sum()
    try:
        bonds = np.atleast_2d(np.loadtxt(name, dtype=int))
    except UserWarning:
        bonds = np.zeros([0,2], int)
    #select bonds that are inside
    bonds = bonds[inside[bonds].min(-1)]
    np.savetxt(
        os.path.join(newpath, os.path.split(name)[-1]),
        bonds,
        fmt='%d'
        )
np.savetxt(
    os.path.join(newpath, x.trajfile[:-4]+'nb'),
    nb, fmt='%d'
    )
x = xp.Experiment(os.path.join(newpath, x.trajfile))
#x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
maxlength = 0
pro = ProgressBar(x.size)

#compute and save the length of broken bonds between consecutive time steps
for t, name in x.enum(ext='bonds'):
    try:
        bonds1 = np.atleast_2d(np.loadtxt(name, dtype=int))
    except UserWarning:
        bonds1 = np.zeros([0,2], int)
    p2tr1 = np.loadtxt(x.get_format_string(ext='p2tr')%t, dtype=int)
    if t>0:
        lengths = broken_bonds_lenghts(bonds0, bonds1, p2tr0, p2tr1)
        np.savetxt(x.get_format_string('_broken', ext='length')%t, lengths, fmt='%d')
        if len(lengths)>0:
            maxlength = max([maxlength, lengths.max()])
    bonds0 = bonds1
    p2tr0 = p2tr1
    pro.animate(t)

#make the time dependent histogram of broken lengths
histlength = np.zeros(maxlength+1, int)
for t, name in x.enum('_broken', ext='length'):
    if t==0: continue
    try:
        lengths = np.loadtxt(name, int)
    except UserWarning:
        continue
    histlength += np.histogram(
        lengths, 
        bins=np.arange(-1, maxlength+1)
        )[0]

#save the time dependent histogram of broken lengths
np.savetxt(
    os.path.join(x.path, 'broken_length.hist'),
    np.column_stack((np.arange(-1, maxlength), histlength)),
    fmt='%d')
x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/2_ageing/1415_ageing.traj')
LL = 128
newpath = os.path.join(x.path, '%d'%LL)
os.mkdir(newpath)
#cheap copy of file.traj
os.symlink(os.path.join(x.path, x.trajfile), os.path.join(newpath, x.trajfile))
#cheap copy of file_t???.p2tr
for t, name in x.enum(ext='p2tr'):
    os.symlink(name, os.path.join(newpath, os.path.split(name)[-1]))

c = 0.5*(m+M)
nb = np.zeros(x.size, int)
#Remove bonds when on of the members are outside
for t, name in x.enum(ext='bonds'):
    pos = np.loadtxt(x.get_format_string()%t, skiprows=2)
    inside = np.min(np.abs(pos-c)<LL/2, -1)
    nb[t] = inside.sum()
    try:
        bonds = np.atleast_2d(np.loadtxt(name, dtype=int))
    except UserWarning:
        bonds = np.zeros([0,2], int)
    #select bonds that are inside
    bonds = bonds[inside[bonds].min(-1)]
    np.savetxt(
        os.path.join(newpath, os.path.split(name)[-1]),
        bonds,
        fmt='%d'
        )
np.savetxt(
    os.path.join(newpath, x.trajfile[:-4]+'nb'),
    nb, fmt='%d'
    )
x = xp.Experiment(os.path.join(newpath, x.trajfile))
#x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
maxlength = 0
pro = ProgressBar(x.size)

#compute and save the length of broken bonds between consecutive time steps
for t, name in x.enum(ext='bonds'):
    try:
        bonds1 = np.atleast_2d(np.loadtxt(name, dtype=int))
    except UserWarning:
        bonds1 = np.zeros([0,2], int)
    p2tr1 = np.loadtxt(x.get_format_string(ext='p2tr')%t, dtype=int)
    if t>0:
        lengths = broken_bonds_lenghts(bonds0, bonds1, p2tr0, p2tr1)
        np.savetxt(x.get_format_string('_broken', ext='length')%t, lengths, fmt='%d')
        if len(lengths)>0:
            maxlength = max([maxlength, lengths.max()])
    bonds0 = bonds1
    p2tr0 = p2tr1
    pro.animate(t)

#make the time dependent histogram of broken lengths
histlength = np.zeros(maxlength+1, int)
for t, name in x.enum('_broken', ext='length'):
    if t==0: continue
    try:
        lengths = np.loadtxt(name, int)
    except UserWarning:
        continue
    histlength += np.histogram(
        lengths, 
        bins=np.arange(-1, maxlength+1)
        )[0]

#save the time dependent histogram of broken lengths
np.savetxt(
    os.path.join(x.path, 'broken_length.hist'),
    np.column_stack((np.arange(-1, maxlength), histlength)),
    fmt='%d')
L = 5
for LL  in ['', '256/', '128/']:
    rates = []
    for head,t0 in zip(['1_percolation/%s1340_percolation', '2_ageing/%s1415_ageing'], [68, 0]):
        x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/%s.traj'%(head%LL))
        lengths = [
            np.loadtxt(name, int) 
            for t, name in x.enum('_broken', 'length') 
            if t>0]
        maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
        his = np.array([
            np.histogram(l, np.arange(-1, maxlength+1))[0]
            for l in lengths])
        rates.append(np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )))
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(li[:-1], [r[i:j].mean() for i,j in zip(li, li[1:])])[0]

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 10*np.arange(1,1000)**(-5.), 'k:')
xscale('log')
yscale('log')
ylim(1e-6,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
L = 5
for LL  in ['', '256/', '128/']:
    rates = []
    for head,t0 in zip(['1_percolation/%s1340_percolation', '2_ageing/%s1415_ageing'], [68, 0]):
        x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/%s.traj'%(head%LL))
        lengths = [
            np.loadtxt(name, int) 
            for t, name in x.enum('_broken', 'length') 
            if t>0]
        maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
        his = np.array([
            np.histogram(l, np.arange(-1, maxlength+1))[0]
            for l in lengths])
        rates.append(np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )))
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])])[0]

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 10*np.arange(1,1000)**(-5.), 'k:')
xscale('log')
yscale('log')
ylim(1e-6,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
M,m
M-m
L = 5
for LL, lab  in zip(['', '256/', '128/'], ['500x500x252', '256x256x252', '128x128x128']):
    rates = []
    for head,t0 in zip(['1_percolation/%s1340_percolation', '2_ageing/%s1415_ageing'], [68, 0]):
        x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/%s.traj'%(head%LL))
        lengths = [
            np.loadtxt(name, int) 
            for t, name in x.enum('_broken', 'length') 
            if t>0]
        maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
        his = np.array([
            np.histogram(l, np.arange(-1, maxlength+1))[0]
            for l in lengths])
        rates.append(np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )))
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], label=lab)[0]

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 10*np.arange(1,1000)**(-5.), 'k:')
xscale('log')
yscale('log')
ylim(1e-6,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
legend('loc=lower right')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/breaking_per_particle_finite_size.pdf')
L = 5
for LL, lab  in zip(['', '256/', '128/'], ['500x500x252', '256x256x252', '128x128x128']):
    rates = []
    for head,t0 in zip(['1_percolation/%s1340_percolation', '2_ageing/%s1415_ageing'], [68, 0]):
        x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/%s.traj'%(head%LL))
        lengths = [
            np.loadtxt(name, int) 
            for t, name in x.enum('_broken', 'length') 
            if t>0]
        maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
        his = np.array([
            np.histogram(l, np.arange(-1, maxlength+1))[0]
            for l in lengths])
        rates.append(np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )))
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], label=lab)[0]

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 10*np.arange(1,1000)**(-5.), 'k:')
xscale('log')
yscale('log')
ylim(1e-6,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
legend(loc='lower right')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/breaking_per_particle_finite_size.pdf')
for L in [5,13, 23]:
    for LL, lab  in zip(['', '256/', '128/'], ['500x500x252', '256x256x252', '128x128x128']):
        rates = []
        for head,t0 in zip(['1_percolation/%s1340_percolation', '2_ageing/%s1415_ageing'], [68, 0]):
            x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/%s.traj'%(head%LL))
            lengths = [
                np.loadtxt(name, int) 
                for t, name in x.enum('_broken', 'length') 
                if t>0]
            maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
            his = np.array([
                np.histogram(l, np.arange(-1, maxlength+1))[0]
                for l in lengths])
            rates.append(np.column_stack((
                np.arange(len(his)-t0),
                his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
                )))
        r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
        #10 points per decade
        li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
        li = np.concatenate(([0],li))
        l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], label=lab)[0]

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 10*np.arange(1,1000)**(-5.), 'k:')
xscale('log')
yscale('log')
ylim(1e-6,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
legend(loc='lower right')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/breaking_per_particle_finite_size_Lthr.pdf')
for L in [5,13, 23]:
    for LL, lab  in zip(['', '256/', '128/'], ['500x500x252', '256x256x252', '128x128x128']):
        rates = []
        for head,t0 in zip(['1_percolation/%s1340_percolation', '2_ageing/%s1415_ageing'], [68, 0]):
            x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/%s.traj'%(head%LL))
            lengths = [
                np.loadtxt(name, int) 
                for t, name in x.enum('_broken', 'length') 
                if t>0]
            maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
            his = np.array([
                np.histogram(l, np.arange(-1, maxlength+1))[0]
                for l in lengths])
            rates.append(np.column_stack((
                np.arange(len(his)-t0),
                his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
                )))
        r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
        #10 points per decade
        li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
        li = np.concatenate(([0],li))
        l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], label=lab)[0]

plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
plot(np.arange(1,1000), 10*np.arange(1,1000)**(-5.), 'k:')
xscale('log')
yscale('log')
ylim(1e-6,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
#legend(loc='lower right')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/breaking_per_particle_finite_size_Lthr.pdf')
x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
lengths = [
    np.loadtxt(name, int) 
    for t, name in x.enum('_broken', 'length') 
    if t>0]
maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
his = np.array([
    np.histogram(l, np.arange(-1, maxlength+1))[0]
    for l in lengths])
plot(np.arange(-1, maxlength), his[69:79].mean(0))
plot(np.arange(-1, maxlength), his[79:109].mean(0))
plot(np.arange(-1, maxlength), his[109:].mean(0))
x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
lengths = [
    np.loadtxt(name, int) 
    for t, name in x.enum('_broken', 'length') 
    if t>0]
maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
his = np.array([
    np.histogram(l, np.arange(-1, maxlength+1))[0]
    for l in lengths])
plot(np.arange(-1, maxlength), his[69:79].mean(0))
plot(np.arange(-1, maxlength), his[79:109].mean(0))
plot(np.arange(-1, maxlength), his[109:].mean(0))
yscale('log')
x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
lengths = [
    np.loadtxt(name, int) 
    for t, name in x.enum('_broken', 'length') 
    if t>0]
maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
his = np.array([
    np.histogram(l, np.arange(-1, maxlength+1))[0]
    for l in lengths])
plot(np.arange(maxlength), his[69:79,2:].mean(0))
plot(np.arange(maxlength), his[79:109,2:].mean(0))
plot(np.arange(maxlength), his[109:,2:].mean(0))
yscale('log')
x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
lengths = [
    np.loadtxt(name, int) 
    for t, name in x.enum('_broken', 'length') 
    if t>0]
maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
his = np.array([
    np.histogram(l, np.arange(-1, maxlength+1))[0]
    for l in lengths])
plot(np.arange(maxlength), his[69:79,1:].mean(0))
plot(np.arange(maxlength), his[79:109,1:].mean(0))
plot(np.arange(maxlength), his[109:,1:].mean(0))
yscale('log')
x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
lengths = [
    np.loadtxt(name, int) 
    for t, name in x.enum('_broken', 'length') 
    if t>0]
maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
his = np.array([
    np.histogram(l, np.arange(-1, maxlength+1))[0]
    for l in lengths])
plot(np.arange(maxlength-1), his[69:79,2:].mean(0))
plot(np.arange(maxlength-1), his[79:109,2:].mean(0))
plot(np.arange(maxlength-1), his[109:,2:].mean(0))
yscale('log')
x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
lengths = [
    np.loadtxt(name, int) 
    for t, name in x.enum('_broken', 'length') 
    if t>0]
maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
his = np.array([
    np.histogram(l, np.arange(-1, maxlength+1))[0]
    for l in lengths])
plot(np.arange(maxlength-1), his[69:79,2:].mean(0))
plot(np.arange(maxlength-1), his[79:109,2:].mean(0))
plot(np.arange(maxlength-1), his[109:,2:].mean(0))
xscale('log');yscale('log')
x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
lengths = [
    np.loadtxt(name, int) 
    for t, name in x.enum('_broken', 'length') 
    if t>0]
maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
his = np.array([
    np.histogram(l, np.arange(-1, maxlength+1))[0]
    for l in lengths])
plot(np.arange(maxlength-1), his[69:79,2:].mean(0))
plot(np.arange(maxlength-1), his[79:109,2:].mean(0))
plot(np.arange(maxlength-1), his[109:,2:].mean(0))
xscale('log');yscale('log')
figure(figsize(24,6))
for i, (LL, lab) in enumerate(zip(['', '256/', '128/'], ['500x500x252', '256x256x252', '128x128x128'])):
    subfigure(3,1,1+i)
    x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/%s1340_percolation.traj'%LL)
    lengths = [
        np.loadtxt(name, int) 
        for t, name in x.enum('_broken', 'length') 
        if t>0]
    maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
    his = np.array([
        np.histogram(l, np.arange(-1, maxlength+1))[0]
        for l in lengths])
    plot(np.arange(maxlength-1), his[69:79,2:].mean(0))
    plot(np.arange(maxlength-1), his[79:109,2:].mean(0))
    plot(np.arange(maxlength-1), his[109:,2:].mean(0))
    xscale('log');yscale('log')
figure(figsize(24,6))
for i, (LL, lab) in enumerate(zip(['', '256/', '128/'], ['500x500x252', '256x256x252', '128x128x128'])):
    subplot(3,1,1+i)
    x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/%s1340_percolation.traj'%LL)
    lengths = [
        np.loadtxt(name, int) 
        for t, name in x.enum('_broken', 'length') 
        if t>0]
    maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
    his = np.array([
        np.histogram(l, np.arange(-1, maxlength+1))[0]
        for l in lengths])
    plot(np.arange(maxlength-1), his[69:79,2:].mean(0))
    plot(np.arange(maxlength-1), his[79:109,2:].mean(0))
    plot(np.arange(maxlength-1), his[109:,2:].mean(0))
    xscale('log');yscale('log')
figure(figsize(24,6))
for i, (LL, lab) in enumerate(zip(['', '256/', '128/'], ['500x500x252', '256x256x252', '128x128x128'])):
    subplot(1,3,1+i)
    x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/%s1340_percolation.traj'%LL)
    lengths = [
        np.loadtxt(name, int) 
        for t, name in x.enum('_broken', 'length') 
        if t>0]
    maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
    his = np.array([
        np.histogram(l, np.arange(-1, maxlength+1))[0]
        for l in lengths])
    plot(np.arange(maxlength-1), his[69:79,2:].mean(0))
    plot(np.arange(maxlength-1), his[79:109,2:].mean(0))
    plot(np.arange(maxlength-1), his[109:,2:].mean(0))
    xscale('log');yscale('log')
figure(figsize(24,6))
ax1 = subplot(1,3,1)
for i, (LL, lab) in enumerate(zip(['', '256/', '128/'], ['500x500x252', '256x256x252', '128x128x128'])):
    if (i>0)
        subplot(1,3,1+i, sharey=ax1)
    x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/%s1340_percolation.traj'%LL)
    lengths = [
        np.loadtxt(name, int) 
        for t, name in x.enum('_broken', 'length') 
        if t>0]
    maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
    his = np.array([
        np.histogram(l, np.arange(-1, maxlength+1))[0]
        for l in lengths])
    plot(np.arange(maxlength-1), his[69:79,2:].mean(0))
    plot(np.arange(maxlength-1), his[79:109,2:].mean(0))
    plot(np.arange(maxlength-1), his[109:,2:].mean(0))
    xscale('log');yscale('log')
figure(figsize(24,6))
ax1 = subplot(1,3,1)
for i, (LL, lab) in enumerate(zip(['', '256/', '128/'], ['500x500x252', '256x256x252', '128x128x128'])):
    if (i>0):
        subplot(1,3,1+i, sharey=ax1)
    x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/%s1340_percolation.traj'%LL)
    lengths = [
        np.loadtxt(name, int) 
        for t, name in x.enum('_broken', 'length') 
        if t>0]
    maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
    his = np.array([
        np.histogram(l, np.arange(-1, maxlength+1))[0]
        for l in lengths])
    plot(np.arange(maxlength-1), his[69:79,2:].mean(0))
    plot(np.arange(maxlength-1), his[79:109,2:].mean(0))
    plot(np.arange(maxlength-1), his[109:,2:].mean(0))
    xscale('log');yscale('log')
figure(figsize(24,6))
ax1 = subplot(1,3,1)
ylabel('pdf')
for i, (LL, lab) in enumerate(zip(['', '256/', '128/'], ['500x500x252', '256x256x252', '128x128x128'])):
    if (i>0):
        subplot(1,3,1+i, sharey=ax1)
    x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/%s1340_percolation.traj'%LL)
    lengths = [
        np.loadtxt(name, int) 
        for t, name in x.enum('_broken', 'length') 
        if t>0]
    maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
    his = np.array([
        np.histogram(l, np.arange(-1, maxlength+1))[0]
        for l in lengths])
    plot(np.arange(maxlength-1), his[69:79,2:].mean(0))
    plot(np.arange(maxlength-1), his[79:109,2:].mean(0))
    plot(np.arange(maxlength-1), his[109:,2:].mean(0))
    xscale('log');yscale('log')
    xlabel('L')
x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/2_ageing/1415_ageing.traj')
bonds0 = np.atleast_2d(np.loadtxt(x.get_format_string(ext='bonds')%0, dtype=int))
bonds1 = np.atleast_2d(np.loadtxt(x.get_format_string(ext='bonds')%1, dtype=int))
p2tr0 = np.loadtxt(x.get_format_string(ext='p2tr')%0, dtype=int)
p2tr1 = np.loadtxt(x.get_format_string(ext='p2tr')%1, dtype=int)
#bonds (between trajectories) existing at t but no more at t+dt
# = broken bonds + lost trajectories
trbonds = set([(a,b) for a,b in np.sort(p2tr0[bonds0], 1)]) - set([(a,b) for a,b in np.sort(p2tr1[bonds1], axis=1)])
#graph of the bonds between trajectories at t+dt
g = nx.Graph()
g.add_nodes_from(p2tr1)
g.add_edges_from(p2tr1[bonds1])
get_ipython().magic(u'pinfo g.add_nodes_from')
bonds0 = np.atleast_2d(np.loadtxt(x.get_format_string(ext='bonds')%0, dtype=int))
bonds1 = np.atleast_2d(np.loadtxt(x.get_format_string(ext='bonds')%1, dtype=int))
p2tr0 = np.loadtxt(x.get_format_string(ext='p2tr')%0, dtype=int)
p2tr1 = np.loadtxt(x.get_format_string(ext='p2tr')%1, dtype=int)
#bonds (between trajectories) existing at t but no more at t+dt
# = broken bonds + lost trajectories
trbonds = set([(a,b) for a,b in np.sort(p2tr0[bonds0], 1)]) - set([(a,b) for a,b in np.sort(p2tr1[bonds1], axis=1)])
#graph of the bonds between trajectories at t+dt
g = nx.Graph()
g.add_nodes_from([(tr, {'p':p}) for p, tr in enumerate(p2tr1)])
g.add_edges_from(p2tr1[bonds1])
trbonds[0]
np.array(trbonds).shape
np.array([[a,b] for a,b in trbonds]).shape
atrbonds = np.array([[a,b] for a,b in trbonds])
nx.shortest_path(g, atrbonds[0,0], atrbonds[0,1])
_[0]
_[0]['p']
atrbonds = np.array([[a,b] for a,b in trbonds])
tr2p1 = dict((tr,p) for p, tr in enumerate(p2tr1))
tr2p1[14761]
get_ipython().magic(u'timeit nx.shortest_path(g, atrbonds[0,0], atrbonds[0,1])')
for a, b in trbonds:
    nx.shortest_path(g, a, b)
for a, b in trbonds:
    if not g.has_node(a) or not g.has_node(b): continue
    nx.shortest_path(g, a, b)
def shortest_path(g, trs):
    """From a graph g, whoes nodes are trajectory indices and edges corresponds to connectivity at time t1,
    and couples of nodes trs, first at time t0 second at time t1,
    return the on-graph shortest path between each couple in term of trajectory indices. 
    
    trs is a iterable of couples of integers."""
    paths = []
    for i, (a,b) in enumerate(trs):
        if not g.has_node(a) or not g.has_node(b): continue
        try:
            paths.append(nx.shortest_path(g, a,b))
        except nx.NetworkXNoPath:
            pass
    return paths
def broken_bonds_path(bonds0, bonds1, p2tr0, p2tr1):
    """Paths on graph at t1 between particles involved in a bond at t0 and no more at t1.
    
    bonds0, bonds1 are respectively the bonds at t0 and 1 in terms of position
    p2tr0, p2tr1 are respectively the position to trajectory relationship at t0 and t1
    
    Returns a list of paths in terms of position indices in t1"""
    #bonds (between trajectories) existing at t but no more at t+dt
    # = broken bonds + lost trajectories
    trbonds = set([(a,b) for a,b in np.sort(p2tr0[bonds0], 1)]) - set([(a,b) for a,b in np.sort(p2tr1[bonds1], axis=1)])

    #graph of the bonds between trajectories at t+dt
    g = nx.Graph()
    g.add_nodes_from(p2tr1)
    g.add_edges_from(p2tr1[bonds1])
    tr2p1 = dict((tr,p) for p, tr in enumerate(p2tr1))
    
    trpaths = shortest_path(g, trbonds)
    ppaths = [
        [tr2p1[tr] for tr in path]
        for path in trpaths
        ]

    return ppaths
ppaths = broken_bonds_path(bonds0, bonds1, p2tr0, p2tr1)
len(ppaths)
def shortest_path(g, trs):
    """From a graph g, whoes nodes are trajectory indices and edges corresponds to connectivity at time t1,
    and couples of nodes trs, first at time t0 second at time t1,
    yield the on-graph shortest path between each couple in term of trajectory indices. 
    
    trs is a iterable of couples of integers."""
    for i, (a,b) in enumerate(trs):
        if not g.has_node(a) or not g.has_node(b): continue
        try:
            yield nx.shortest_path(g, a,b)
        except nx.NetworkXNoPath:
            pass
def broken_bonds_path(bonds0, bonds1, p2tr0, p2tr1):
    """Paths on graph at t1 between particles involved in a bond at t0 and no more at t1.
    
    bonds0, bonds1 are respectively the bonds at t0 and 1 in terms of position
    p2tr0, p2tr1 are respectively the position to trajectory relationship at t0 and t1
    
    generate paths in terms of position indices in t1"""
    #bonds (between trajectories) existing at t but no more at t+dt
    # = broken bonds + lost trajectories
    trbonds = set([(a,b) for a,b in np.sort(p2tr0[bonds0], 1)]) - set([(a,b) for a,b in np.sort(p2tr1[bonds1], axis=1)])

    #graph of the bonds between trajectories at t+dt
    g = nx.Graph()
    g.add_nodes_from(p2tr1)
    g.add_edges_from(p2tr1[bonds1])
    tr2p1 = dict((tr,p) for p, tr in enumerate(p2tr1))
    
    for trpath in shortest_path(g, trbonds):
        yield [tr2p1[tr] for tr in trpath]
ppaths = list(broken_bonds_path(bonds0, bonds1, p2tr0, p2tr1))
len(ppaths)
get_ipython().magic(u'pinfo %mlab')
g.number_of_edges()
get_ipython().magic(u'pinfo g.number_of_edges')
nbpass = np.zeros(g.number_of_edges(), int)
for trpath in shortest_path(g, trbonds):
    nbpass[trpath] += 1
h, bins = np.histogram(nbpass)
plot(bins[:-1], h)
nbpass = np.zeros(g.number_of_edges(), int)
for trpath in shortest_path(g, trbonds):
    nbpass[trpath] += 1
h, bins = np.histogram(nbpass)
step(bins[:-1], h)
nbpass = np.zeros(g.number_of_edges(), int)
for trpath in shortest_path(g, trbonds):
    nbpass[trpath] += 1
h, bins = np.histogram(nbpass)
step(bins[:-1], h)
yscale('log')
nbpass = np.zeros(g.number_of_nodes(), int)
for trpath in shortest_path(g, trbonds):
    nbpass[trpath] += 1
h, bins = np.histogram(nbpass)
step(bins[:-1], h)
yscale('log')
nbpass = np.zeros(g.number_of_nodes(), int)
for trpath in shortest_path(g, trbonds):
    nbpass[trpath] += 1
h, bins = np.histogram(nbpass, bins=np.arange(5))
step(bins[:-1], h)
yscale('log')
nbpass = np.zeros(g.number_of_nodes(), int)
for trpath in shortest_path(g, trbonds):
    nbpass[trpath] += 1
h, bins = np.histogram(nbpass, bins=np.arange(10)
step(bins[:-1], h)
yscale('log')
nbpass = np.zeros(g.number_of_nodes(), int)
for trpath in shortest_path(g, trbonds):
    nbpass[trpath] += 1
h, bins = np.histogram(nbpass, bins=np.arange(10))
step(bins[:-1], h)
yscale('log')
nbpass = np.zeros(len(p2tr1), int)
for trpath in ppaths:
    nbpass[trpath] += 1
h, bins = np.histogram(nbpass, bins=np.arange(10))
step(bins[:-1], h)
yscale('log')
np.where(nbpass>2)
np.where(nbpass>3)
for ppath in ppaths:
    if 87 in path: print path
for path in ppaths:
    if 87 in path: print path
np.where(nbpass>2)
for path in ppaths:
    if 19 in path: print path
for path in ppaths:
    if 26 in path: print path
for path in ppaths:
    if 118 in path: print path
for path in ppaths:
    if 19 in path: print path
for path in ppaths:
    if 426 in path: print path
for path in ppaths:
    if 677 in path: print path
for path in ppaths:
    if 3223 in path: print path
for path in ppaths:
    if 3282 in path: print path
for path in ppaths:
    if 19575 in path: print path
for path in ppaths:
    if 25957 in path: print path
from colloids import vtk
v = vtk.Polydata()
pos1 = np.loadtxt(x.get_format_string()%1, skiprows=2)
v.points = pos1
v.bonds = np.array([[a,b] for a,b in set((i,j) for i,j in np.sort([(u,v) for path in ppaths for u,v in zip(path, path[1:])], axis=-1))])
v = vtk.Polydata()
pos1 = np.loadtxt(x.get_format_string()%1, skiprows=2)
v.points = pos1
v.bonds = np.array([[a,b] for a,b in set((i,j) for i,j in np.sort([(k,l) for path in ppaths for k,l in zip(path, path[1:])], axis=-1))])
v.bonds.shape
v.bonds.dtype
v.save(x.get_format_string('_brokenpaths', 'vtk')%1)
x.get_format_string('_brokenpaths', 'vtk')%1
v = vtk.Polydata()
pos1 = np.loadtxt(x.get_format_string()%1, skiprows=2)
v.points = pos1
v.bonds = np.array([[a,b] for a,b in set((i,j) for i,j in np.sort([(k,l) for path in ppaths for k,l in zip(path, path[1:])], axis=-1))])
lengths = np.zeros(len(v.bonds), int)
for path in ppaths:
    lengths[path] = np.maximum(lengths[path], len(path))
v.bondsScalars = ['length', lengths]
nbpass = np.zeros(len(pos1), int)
for path in ppaths:
    nbpass[path] += 1
v.scalars = ['nbpass', nbpass]
v.save(x.get_format_string('_brokenpaths', 'vtk')%1)
v = vtk.Polydata()
pos1 = np.loadtxt(x.get_format_string()%1, skiprows=2)
v.points = pos1
v.bonds = np.array([[a,b] for a,b in set((i,j) for i,j in np.sort([(k,l) for path in ppaths for k,l in zip(path, path[1:])], axis=-1))])
lengths = np.zeros(len(v.bonds), int)
for path in ppaths:
    for k,l in zip(path, path[1:]):
        i,j = np.sort([k,l])
        b = np.where(v.bonds == [i,j])[0]
        lengths[b] = max([lengths[b], len(path)])
v.bondsScalars = ['length', lengths]
nbpass = np.zeros(len(pos1), int)
for path in ppaths:
    nbpass[path] += 1
v.scalars = ['nbpass', nbpass]
v.save(x.get_format_string('_brokenpaths', 'vtk')%1)
v = vtk.Polydata()
pos1 = np.loadtxt(x.get_format_string()%1, skiprows=2)
v.points = pos1
v.bonds = np.array([[a,b] for a,b in set((i,j) for i,j in np.sort([(k,l) for path in ppaths for k,l in zip(path, path[1:])], axis=-1))])
lengths = np.zeros(len(v.bonds), int)
for path in ppaths:
    for k,l in zip(path, path[1:]):
        i,j = np.sort([k,l])
        b = np.where(v.bonds[:,0] == i,j])[0][0]
        lengths[b] = max([lengths[b], len(path)])
v.bondsScalars = ['length', lengths]
nbpass = np.zeros(len(pos1), int)
for path in ppaths:
    nbpass[path] += 1
v.scalars = ['nbpass', nbpass]
v.save(x.get_format_string('_brokenpaths', 'vtk')%1)
v = vtk.Polydata()
pos1 = np.loadtxt(x.get_format_string()%1, skiprows=2)
v.points = pos1
v.bonds = np.array([[a,b] for a,b in set((i,j) for i,j in np.sort([(k,l) for path in ppaths for k,l in zip(path, path[1:])], axis=-1))])
lengths = np.zeros(len(v.bonds), int)
for path in ppaths:
    for k,l in zip(path, path[1:]):
        i,j = np.sort([k,l])
        b = np.where(v.bonds[:,0] == [i,j])[0][0]
        lengths[b] = max([lengths[b], len(path)])
v.bondsScalars = ['length', lengths]
nbpass = np.zeros(len(pos1), int)
for path in ppaths:
    nbpass[path] += 1
v.scalars = ['nbpass', nbpass]
v.save(x.get_format_string('_brokenpaths', 'vtk')%1)
get_ipython().magic(u'debug ')
v = vtk.Polydata()
pos1 = np.loadtxt(x.get_format_string()%1, skiprows=2)
v.points = pos1
v.bonds = np.array([[a,b] for a,b in set((i,j) for i,j in np.sort([(k,l) for path in ppaths for k,l in zip(path, path[1:])], axis=-1))])
lengths = np.zeros(len(v.bonds), int)
for path in ppaths:
    for k,l in zip(path, path[1:]):
        i,j = np.sort([k,l])
        b = np.where(v.bonds == [i,j])[0][0]
        lengths[b] = max([lengths[b], len(path)])
v.bondsScalars = ['length', lengths]
nbpass = np.zeros(len(pos1), int)
for path in ppaths:
    nbpass[path] += 1
v.scalars = ['nbpass', nbpass]
v.save(x.get_format_string('_brokenpaths', 'vtk')%1)
v = vtk.Polydata()
pos1 = np.loadtxt(x.get_format_string()%1, skiprows=2)
v.points = pos1
v.bonds = np.array([[a,b] for a,b in set((i,j) for i,j in np.sort([(k,l) for path in ppaths for k,l in zip(path, path[1:])], axis=-1))])
lengths = np.zeros(len(v.bonds), int)
for path in ppaths:
    for k,l in zip(path, path[1:]):
        i,j = np.sort([k,l])
        b = np.where(v.bonds == [i,j])[0][0]
        lengths[b] = max([lengths[b], len(path)])
v.bondsScalars = [('length', lengths)]
nbpass = np.zeros(len(pos1), int)
for path in ppaths:
    nbpass[path] += 1
v.scalars = [('nbpass', nbpass)]
v.save(x.get_format_string('_brokenpaths', 'vtk')%1)
x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
maxlength = 0
pro = ProgressBar(x.size)

#visualise broken bonds between consecutive time steps
for t, name in x.enum(ext='bonds'):
    try:
        bonds1 = np.atleast_2d(np.loadtxt(name, dtype=int))
    except UserWarning:
        bonds1 = np.zeros([0,2], int)
    p2tr1 = np.loadtxt(x.get_format_string(ext='p2tr')%t, dtype=int)
    if t>0:
        ppaths = list(broken_bonds_path(bonds0, bonds1, p2tr0, p2tr1))
        #export all positions
        v = vtk.Polydata()
        pos1 = np.loadtxt(x.get_format_string()%1, skiprows=2)
        v.points = pos1
        #export broken paths
        v.bonds = np.array([[a,b] for a,b in set((i,j) for i,j in np.sort([(k,l) for path in ppaths for k,l in zip(path, path[1:])], axis=-1))])
        #label broken path by their lengths
        lengths = np.zeros(len(v.bonds), int)
        for path in ppaths:
            for k,l in zip(path, path[1:]):
                i,j = np.sort([k,l])
                b = np.where(v.bonds == [i,j])[0][0]
                lengths[b] = max([lengths[b], len(path)])
        v.bondsScalars = [('length', lengths)]
        #label particles by how many broken paths they participate in
        nbpass = np.zeros(len(pos1), int)
        for path in ppaths:
            nbpass[path] += 1
        v.scalars = [('nbpass', nbpass)]
        #to disk
        v.save(x.get_format_string('_brokenpaths', 'vtk')%1)
    bonds0 = bonds1
    p2tr0 = p2tr1
    pro.animate(t)
x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
maxlength = 0
pro = ProgressBar(x.size)

#visualise broken bonds between consecutive time steps
for t, name in x.enum(ext='bonds'):
    try:
        bonds1 = np.atleast_2d(np.loadtxt(name, dtype=int))
    except UserWarning:
        bonds1 = np.zeros([0,2], int)
    p2tr1 = np.loadtxt(x.get_format_string(ext='p2tr')%t, dtype=int)
    if t>0:
        ppaths = list(broken_bonds_path(bonds0, bonds1, p2tr0, p2tr1))
        #export all positions
        v = vtk.Polydata()
        pos1 = np.loadtxt(x.get_format_string()%t, skiprows=2)
        v.points = pos1
        #export broken paths
        v.bonds = np.array([[a,b] for a,b in set((i,j) for i,j in np.sort([(k,l) for path in ppaths for k,l in zip(path, path[1:])], axis=-1))])
        #label broken path by their lengths
        lengths = np.zeros(len(v.bonds), int)
        for path in ppaths:
            for k,l in zip(path, path[1:]):
                i,j = np.sort([k,l])
                b = np.where(v.bonds == [i,j])[0][0]
                lengths[b] = max([lengths[b], len(path)])
        v.bondsScalars = [('length', lengths)]
        #label particles by how many broken paths they participate in
        nbpass = np.zeros(len(pos1), int)
        for path in ppaths:
            nbpass[path] += 1
        v.scalars = [('nbpass', nbpass)]
        #to disk
        v.save(x.get_format_string('_brokenpaths', 'vtk')%t)
    bonds0 = bonds1
    p2tr0 = p2tr1
    pro.animate(t)
#x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
maxlength = 0
pro = ProgressBar(x.size)

#compute and save the length of broken bonds between consecutive time steps
for t, name in x.enum(ext='bonds'):
    try:
        bonds1 = np.atleast_2d(np.loadtxt(name, dtype=int))
    except UserWarning:
        bonds1 = np.zeros([0,2], int)
    p2tr1 = np.loadtxt(x.get_format_string(ext='p2tr')%t, dtype=int)
    if t>0:
        lengths = []
        pos1 = np.loadtxt(x.get_format_string()%t, skiprows=2)
        #bounding box
        m = pos1.min(-1)
        M = pos1.max(-1)
        #filter broken bonds' shortest path
        for path in broken_bonds_path(bonds0, bonds1, p2tr0, p2tr1):
            #positions in the path
            pos = pos1[path]
            #lower boundaries
            #ex-bond members closest to the lower boundaries
            p = pos[[0,-1]][np.argmin(pos[[0,-1]]-m, axis=0)]
            #path extent away from this position
            extent = np.max(pos - p, axis=0)
            #is the path further than its own extent from the boundaries?
            if not min(p - extent < m):
                continue
            #upper boundaries
            #ex-bond members closest to the upper boundaries
            p = pos[[0,-1]][np.argmin(M - pos[[0,-1]], axis=0)]
            #path extent away from this position
            extent = np.max(p - pos, axis=0)
            #is the path further than its own extent from the boundaries?
            if not min(p + extent > M):
                continue
            lengths.append(len(path))
            
        np.savetxt(x.get_format_string('_brokenbulk', ext='length')%t, lengths, fmt='%d')
        if len(lengths)>0:
            maxlength = max([maxlength, max(lengths)])
    bonds0 = bonds1
    p2tr0 = p2tr1
    pro.animate(t)

#make the time dependent histogram of broken lengths
histlength = np.zeros(maxlength+1, int)
for t, name in x.enum('_brokenbulk', ext='length'):
    if t==0: continue
    try:
        lengths = np.loadtxt(name, int)
    except UserWarning:
        continue
    histlength += np.histogram(
        lengths, 
        bins=np.arange(-1, maxlength+1)
        )[0]

#save the time dependent histogram of broken lengths
np.savetxt(
    os.path.join(x.path, 'brokenbulk_length.hist'),
    np.column_stack((np.arange(-1, maxlength), histlength)),
    fmt='%d')
get_ipython().magic(u'debug ')
#x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
maxlength = 0
pro = ProgressBar(x.size)

#compute and save the length of broken bonds between consecutive time steps
for t, name in x.enum(ext='bonds'):
    try:
        bonds1 = np.atleast_2d(np.loadtxt(name, dtype=int))
    except UserWarning:
        bonds1 = np.zeros([0,2], int)
    p2tr1 = np.loadtxt(x.get_format_string(ext='p2tr')%t, dtype=int)
    if t>0:
        lengths = []
        pos1 = np.loadtxt(x.get_format_string()%t, skiprows=2)
        #bounding box
        m = pos1.min(0)
        M = pos1.max(0)
        #filter broken bonds' shortest path
        for path in broken_bonds_path(bonds0, bonds1, p2tr0, p2tr1):
            #positions in the path
            pos = pos1[path]
            #lower boundaries
            #ex-bond members closest to the lower boundaries
            p = pos[[0,-1]][np.argmin(pos[[0,-1]]-m, axis=0)]
            #path extent away from this position
            extent = np.max(pos - p, axis=0)
            #is the path further than its own extent from the boundaries?
            if not min(p - extent < m):
                continue
            #upper boundaries
            #ex-bond members closest to the upper boundaries
            p = pos[[0,-1]][np.argmin(M - pos[[0,-1]], axis=0)]
            #path extent away from this position
            extent = np.max(p - pos, axis=0)
            #is the path further than its own extent from the boundaries?
            if not min(p + extent > M):
                continue
            lengths.append(len(path))
            
        np.savetxt(x.get_format_string('_brokenbulk', ext='length')%t, lengths, fmt='%d')
        if len(lengths)>0:
            maxlength = max([maxlength, max(lengths)])
    bonds0 = bonds1
    p2tr0 = p2tr1
    pro.animate(t)

#make the time dependent histogram of broken lengths
histlength = np.zeros(maxlength+1, int)
for t, name in x.enum('_brokenbulk', ext='length'):
    if t==0: continue
    try:
        lengths = np.loadtxt(name, int)
    except UserWarning:
        continue
    histlength += np.histogram(
        lengths, 
        bins=np.arange(-1, maxlength+1)
        )[0]

#save the time dependent histogram of broken lengths
np.savetxt(
    os.path.join(x.path, 'brokenbulk_length.hist'),
    np.column_stack((np.arange(-1, maxlength), histlength)),
    fmt='%d')
get_ipython().magic(u'debug ')
#x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
maxlength = 0
pro = ProgressBar(x.size)

#compute and save the length of broken bonds between consecutive time steps
for t, name in x.enum(ext='bonds'):
    try:
        bonds1 = np.atleast_2d(np.loadtxt(name, dtype=int))
    except UserWarning:
        bonds1 = np.zeros([0,2], int)
    p2tr1 = np.loadtxt(x.get_format_string(ext='p2tr')%t, dtype=int)
    if t>0:
        lengths = []
        pos1 = np.loadtxt(x.get_format_string()%t, skiprows=2)
        #bounding box
        m = pos1.min(0)
        M = pos1.max(0)
        #filter broken bonds' shortest path
        for path in broken_bonds_path(bonds0, bonds1, p2tr0, p2tr1):
            #positions in the path
            pos = pos1[path]
            #lower boundaries
            #ex-bond members closest to the lower boundaries
            p = m + np.min(pos[[0,-1]]-m, axis=0)
            #path extent away from this position
            extent = np.max(pos - p, axis=0)
            #is the path further than its own extent from the boundaries?
            if not min(p - extent < m):
                continue
            #upper boundaries
            #ex-bond members closest to the upper boundaries
            p = M - np.min(M - pos[[0,-1]], axis=0)
            #path extent away from this position
            extent = np.max(p - pos, axis=0)
            #is the path further than its own extent from the boundaries?
            if not min(p + extent > M):
                continue
            lengths.append(len(path))
            
        np.savetxt(x.get_format_string('_brokenbulk', ext='length')%t, lengths, fmt='%d')
        if len(lengths)>0:
            maxlength = max([maxlength, max(lengths)])
    bonds0 = bonds1
    p2tr0 = p2tr1
    pro.animate(t)

#make the time dependent histogram of broken lengths
histlength = np.zeros(maxlength+1, int)
for t, name in x.enum('_brokenbulk', ext='length'):
    if t==0: continue
    try:
        lengths = np.loadtxt(name, int)
    except UserWarning:
        continue
    histlength += np.histogram(
        lengths, 
        bins=np.arange(-1, maxlength+1)
        )[0]

#save the time dependent histogram of broken lengths
np.savetxt(
    os.path.join(x.path, 'brokenbulk_length.hist'),
    np.column_stack((np.arange(-1, maxlength), histlength)),
    fmt='%d')
p = m + np.min(pos[[0,-1]]-m, axis=0)
p = m + np.min(pos[[0,-1]]-m, axis=0)
print p
extent = np.max(pos - p, axis=0)
print extent
min(p - extent < m)
p = m + np.min(pos[[0,-1]]-m, axis=0)
print p
extent = np.max(pos - p, axis=0)
print extent
max(p - extent < m)
p = M - np.min(M - pos[[0,-1]], axis=0)
print p
extent = np.max(p - pos, axis=0)
print extent
max(p + extent > M)
#x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
maxlength = 0
pro = ProgressBar(x.size)

#compute and save the length of broken bonds between consecutive time steps
for t, name in x.enum(ext='bonds'):
    try:
        bonds1 = np.atleast_2d(np.loadtxt(name, dtype=int))
    except UserWarning:
        bonds1 = np.zeros([0,2], int)
    p2tr1 = np.loadtxt(x.get_format_string(ext='p2tr')%t, dtype=int)
    if t>0:
        lengths = []
        pos1 = np.loadtxt(x.get_format_string()%t, skiprows=2)
        #bounding box
        m = pos1.min(0)
        M = pos1.max(0)
        #filter broken bonds' shortest path
        for path in broken_bonds_path(bonds0, bonds1, p2tr0, p2tr1):
            #positions in the path
            pos = pos1[path]
            #lower boundaries
            #ex-bond members closest to the lower boundaries
            p = m + np.min(pos[[0,-1]]-m, axis=0)
            #path extent away from this position
            extent = np.max(pos - p, axis=0)
            #is the path further than its own extent from the boundaries?
            if max(p - extent < m):
                continue
            #upper boundaries
            #ex-bond members closest to the upper boundaries
            p = M - np.min(M - pos[[0,-1]], axis=0)
            #path extent away from this position
            extent = np.max(p - pos, axis=0)
            #is the path further than its own extent from the boundaries?
            if max(p + extent > M):
                continue
            lengths.append(len(path))
            
        np.savetxt(x.get_format_string('_brokenbulk', ext='length')%t, lengths, fmt='%d')
        if len(lengths)>0:
            maxlength = max([maxlength, max(lengths)])
    bonds0 = bonds1
    p2tr0 = p2tr1
    pro.animate(t)

#make the time dependent histogram of broken lengths
histlength = np.zeros(maxlength+1, int)
for t, name in x.enum('_brokenbulk', ext='length'):
    if t==0: continue
    try:
        lengths = np.loadtxt(name, int)
    except UserWarning:
        continue
    histlength += np.histogram(
        lengths, 
        bins=np.arange(-1, maxlength+1)
        )[0]

#save the time dependent histogram of broken lengths
np.savetxt(
    os.path.join(x.path, 'brokenbulk_length.hist'),
    np.column_stack((np.arange(-1, maxlength), histlength)),
    fmt='%d')
plot(np.arange(len(histlength)-t0), histlength[t0:,4:].sum(-1))
maxlength
histlength.shape
plot(histlength)
lengths = [
    np.loadtxt(name, int) 
    for t, name in x.enum('_brokenbulk', 'length') 
    if t>0]
maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
hist = np.array([
    np.histogram(l, np.arange(maxlength+1))[0]
    for l in lengths])
plot(hist)
lengths = [
    np.loadtxt(name, int) 
    for t, name in x.enum('_brokenbulk', 'length') 
    if t>0]
maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
hist = np.array([
    np.histogram(l, np.arange(maxlength+1))[0]
    for l in lengths])
plot(hist[t0:].mean(0))
lengths = [
    np.loadtxt(name, int) 
    for t, name in x.enum('_brokenbulk', 'length') 
    if t>0]
maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
hist = np.array([
    np.histogram(l, np.arange(maxlength+1))[0]
    for l in lengths])
plot(hist[t0:].mean(0))
xscale('log');yscale('log')
len(hists)
plot(hists[0][t0:].mean(0))
plot(hist[t0:].mean(0))
xscale('log');yscale('log')
plot(np.arange(-1, hists[0].shape[1]-1), hists[0][t0:].mean(0))
plot(hist[t0:].mean(0))
xscale('log');yscale('log')
lengths = [
    np.loadtxt(name, int) 
    for t, name in x.enum('_brokenbulk', 'length') 
    if t>0]
maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
hist = np.array([
    np.histogram(l, np.arange(-1, maxlength+1))[0]
    for l in lengths])
plot(hist[t0:].mean(0))
xscale('log');yscale('log')
plot(np.arange(-1, hists[0].shape[1]-1), hists[0][t0:].mean(0))
plot(np.arange(-1, hist.shape[1]-1), hist[t0:].mean(0))
xscale('log');yscale('log')
hitsts[0][-1,:5]
hists[0][-1,:5]
print hists[0][-1,:5]
print hist[-1,:5]
x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
maxlength = 0
pro = ProgressBar(x.size)

#visualise broken bonds between consecutive time steps
for t, name in x.enum(ext='bonds'):
    try:
        bonds1 = np.atleast_2d(np.loadtxt(name, dtype=int))
    except UserWarning:
        bonds1 = np.zeros([0,2], int)
    p2tr1 = np.loadtxt(x.get_format_string(ext='p2tr')%t, dtype=int)
    if t>0:
        ppaths = list(broken_bonds_path(bonds0, bonds1, p2tr0, p2tr1))
        #export all positions
        v = vtk.Polydata()
        pos1 = np.loadtxt(x.get_format_string()%t, skiprows=2)
        v.points = pos1
        #export broken paths
        v.bonds = np.array([[a,b] for a,b in set((i,j) for i,j in np.sort([(k,l) for path in ppaths for k,l in zip(path, path[1:])], axis=-1))])
        #label broken path by their lengths
        lengths = np.zeros(len(v.bonds), int)
        for path in ppaths:
            for k,l in zip(path, path[1:]):
                i,j = np.sort([k,l])
                b = np.where(v.bonds == [i,j])[0][0]
                lengths[b] = max([lengths[b], len(path)-1])
        v.bondsScalars = [('length', lengths)]
        #label particles by how many broken paths they participate in
        nbpass = np.zeros(len(pos1), int)
        for path in ppaths:
            nbpass[path] += 1
        v.scalars = [('nbpass', nbpass)]
        #to disk
        v.save(x.get_format_string('_brokenpaths', 'vtk')%t)
    bonds0 = bonds1
    p2tr0 = p2tr1
    pro.animate(t)
#x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
maxlength = 0
pro = ProgressBar(x.size)

#compute and save the length of broken bonds between consecutive time steps
for t, name in x.enum(ext='bonds'):
    try:
        bonds1 = np.atleast_2d(np.loadtxt(name, dtype=int))
    except UserWarning:
        bonds1 = np.zeros([0,2], int)
    p2tr1 = np.loadtxt(x.get_format_string(ext='p2tr')%t, dtype=int)
    if t>0:
        lengths = []
        pos1 = np.loadtxt(x.get_format_string()%t, skiprows=2)
        #bounding box
        m = pos1.min(0)
        M = pos1.max(0)
        #filter broken bonds' shortest path
        for path in broken_bonds_path(bonds0, bonds1, p2tr0, p2tr1):
            #positions in the path
            pos = pos1[path]
            #lower boundaries
            #ex-bond members closest to the lower boundaries
            p = m + np.min(pos[[0,-1]]-m, axis=0)
            #path extent away from this position
            extent = np.max(pos - p, axis=0)
            #is the path further than its own extent from the boundaries?
            if max(p - extent < m):
                continue
            #upper boundaries
            #ex-bond members closest to the upper boundaries
            p = M - np.min(M - pos[[0,-1]], axis=0)
            #path extent away from this position
            extent = np.max(p - pos, axis=0)
            #is the path further than its own extent from the boundaries?
            if max(p + extent > M):
                continue
            lengths.append(len(path)-1)
            
        np.savetxt(x.get_format_string('_brokenbulk', ext='length')%t, lengths, fmt='%d')
        if len(lengths)>0:
            maxlength = max([maxlength, max(lengths)])
    bonds0 = bonds1
    p2tr0 = p2tr1
    pro.animate(t)

#make the time dependent histogram of broken lengths
histlength = np.zeros((x.size-1, maxlength+1), int)
for t, name in x.enum('_brokenbulk', ext='length'):
    if t==0: continue
    try:
        lengths = np.loadtxt(name, int)
    except UserWarning:
        continue
    histlength[t-1] = np.histogram(
        lengths, 
        bins=np.arange(-1, maxlength+1)
        )[0]

#save the time dependent histogram of broken lengths
np.savetxt(
    os.path.join(x.path, 'brokenbulk_length.hist'),
    np.column_stack((np.arange(-1, maxlength), histlength)),
    fmt='%d')
lengths = [
    np.loadtxt(name, int) 
    for t, name in x.enum('_brokenbulk', 'length') 
    if t>0]
maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
hist = np.array([
    np.histogram(l, np.arange(-1, maxlength+1))[0]
    for l in lengths])
plot(hist[t0:].mean(0))
xscale('log');yscale('log')
plot(np.arange(-1, hists[0].shape[1]-1), hists[0][t0:].mean(0))
plot(np.arange(-1, hist.shape[1]-1), hist[t0:].mean(0))
xscale('log');yscale('log')
print hists[0][-1,:5]
print hist[-1,:5]
rates = [h[t0:, 5:].sum()*1./x.get_nb()[t0:-1] for h in [hists[0], hist]]
for r in rates:
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(li[:-1], [r[i:j].mean() for i,j in zip(li, li[1:])], label=lab)[0]
xscale('log')
yscale('log')
ylim(1e-6,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
rates = [h[t0:, 5:].sum()*1./x.get_nb()[t0:-1] for h in [hists[0], hist]]
for r in rates:
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(li[:-1], [r[i:j].mean() for i,j in zip(li, li[1:])], label=lab)[0]
xscale('log')
yscale('log')
#ylim(1e-6,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
for na in ['', 'bulk']:    
    lengths = [
        np.loadtxt(name, int) 
        for t, name in x.enum('_broken'+na, 'length') 
        if t>0]
    maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
    his = np.array([
        np.histogram(l, np.arange(-1, maxlength+1))[0]
        for l in lengths])
    r = his[t0:, 5:].sum()*1./x.get_nb()[t0:-1]
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(li[:-1], [r[i:j].mean() for i,j in zip(li, li[1:])], label=lab)[0]
xscale('log')
yscale('log')
#ylim(1e-6,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
L=5
for na in ['', 'bulk']:    
    lengths = [
        np.loadtxt(name, int) 
        for t, name in x.enum('_broken'+na, 'length') 
        if t>0]
    maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
    his = np.array([
        np.histogram(l, np.arange(-1, maxlength+1))[0]
        for l in lengths])
    r = his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(li[:-1], [r[i:j].mean() for i,j in zip(li, li[1:])], label=lab)[0]
xscale('log')
yscale('log')
#ylim(1e-6,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
L=5
t0 = 68
for na in ['', 'bulk']:    
    lengths = [
        np.loadtxt(name, int) 
        for t, name in x.enum('_broken'+na, 'length') 
        if t>0]
    maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
    his = np.array([
        np.histogram(l, np.arange(-1, maxlength+1))[0]
        for l in lengths])
    r = his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(li[:-1], [r[i:j].mean() for i,j in zip(li, li[1:])], label=lab)[0]
xscale('log')
yscale('log')
#ylim(1e-6,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
L=5
t0 = 68
for na in ['', 'bulk']:    
    lengths = [
        np.loadtxt(name, int) 
        for t, name in x.enum('_broken'+na, 'length') 
        if t>0]
    plot([l[l>2].mean() for l in lengths[t0:]])
xscale('log')
xlabel(r'$t/\tau_B$')
ylabel(r'<L>')
t0 = 68
figure(figsize=(16,6))
for na in ['', 'bulk']:    
    lengths = [
        np.loadtxt(name, int) 
        for t, name in x.enum('_broken'+na, 'length') 
        if t>0]
    for i,L in enumerate([2, 10, 20]):
        subplot(3,1,1+i).plot([l[l>L].mean() for l in lengths[t0:]])
xscale('log')
xlabel(r'$t/\tau_B$')
ylabel(r'<L>')
t0 = 68
figure(figsize=(16,6))
for na in ['', 'bulk']:    
    lengths = [
        np.loadtxt(name, int) 
        for t, name in x.enum('_broken'+na, 'length') 
        if t>0]
    for i,L in enumerate([2, 10, 20]):
        subplot(3,1,1+i).plot([l[l>L].mean() for l in lengths[t0:]])
for i in range(3):
    subplot(3,1,1+i)
    xscale('log')
    ylabel(r'<L>')
xlabel(r'$t/\tau_B$')
t0 = 68
figure(figsize=(16,6))
ax1 = subplot(1,3,1)
for na in ['', 'bulk']:    
    lengths = [
        np.loadtxt(name, int) 
        for t, name in x.enum('_broken'+na, 'length') 
        if t>0]
    for i,L in enumerate([2, 10, 20]):
        subplot(3,1,1+i, sharey=ax1).plot([l[l>L].mean() for l in lengths[t0:]])
for i in range(3):
    subplot(1,3,1+i)
    xscale('log')
    ylabel(r'<L>')
xlabel(r'$t/\tau_B$')
t0 = 68
figure(figsize=(16,6))
ax1 = subplot(1,3,1)
for na in ['', 'bulk']:    
    lengths = [
        np.loadtxt(name, int) 
        for t, name in x.enum('_broken'+na, 'length') 
        if t>0]
    for i,L in enumerate([2, 10, 20]):
        subplot(3,1,1+i, sharey=ax1).plot([l[l>L].mean() for l in lengths[t0:]])
for i in range(3):
    subplot(1,3,1+i, sharey=ax1)
    xscale('log')
    ylabel(r'<L>')
xlabel(r'$t/\tau_B$')
t0 = 68
figure(figsize=(16,6))
ax1 = subplot(1,3,1)
xscale('log')
ax2 = subplot(1,3,2, sharey=ax1)
xscale('log')
ax3 = subplot(1,3,3, sharey=ax1)
xscale('log')
for na in ['', 'bulk']:    
    lengths = [
        np.loadtxt(name, int) 
        for t, name in x.enum('_broken'+na, 'length') 
        if t>0]
    for ax,L in zip([ax1, ax2, ax3],[2, 10, 20]):
        ax.plot([l[l>L].mean() for l in lengths[t0:]])
        xlabel(r'$t/\tau_B$')
t0 = 68
figure(figsize=(16,6))
ax1 = subplot(1,3,1)
ylabel(r'<L>')
xscale('log')
ax2 = subplot(1,3,2, sharey=ax1)
xscale('log')
ax3 = subplot(1,3,3, sharey=ax1)
xscale('log')
for na in ['', 'bulk']:    
    lengths = [
        np.loadtxt(name, int) 
        for t, name in x.enum('_broken'+na, 'length') 
        if t>0]
    for ax,L in zip([ax1, ax2, ax3],[2, 10, 20]):
        ax.plot([l[l>L].mean() for l in lengths[t0:]])
        xlabel(r'$t/\tau_B$')
for i, t in enumerate(t0 +2**np.arange(5)):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    h,bins = np.histogram(lengths, bins = np.arange(-1, max(lengths)+1))
    plot(bins[:-1], h, color=cm.autumn(i/5.))
for i, t in enumerate(t0 +2**np.arange(5)):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    h,bins = np.histogram(lengths, bins = np.arange(-1, max(lengths)+1))
    plot(bins[:-1], h, color=cm.autumn(i/5.))
xscale('log'); yscale('log')
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0), 10*np.log10(x.size-t0)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    h,bins = np.histogram(lengths, bins = np.arange(-1, max(lengths)+1))
    plot(bins[:-1], h, color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    h,bins = np.histogram(lengths, bins = np.arange(-1, max(lengths)+1))
    plot(bins[:-1], h, color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    bonds = np.loadtxt(x.get_format_string(ext='bonds')%t, dtype=int)
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    h,bins = np.histogram(lengths, bins = np.arange(-1, max(lengths)+1))
    plot(bins[:-1], h*1./len(bonds), color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L'); ylabel('pdf/#bonds')
f = gcf()
f.aname
f.name
f.get_label()
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    bonds = np.loadtxt(x.get_format_string(ext='bonds')%t, dtype=int)
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    h,bins = np.histogram(lengths, bins = np.arange(-1, max(lengths)+1))
    plot(bins[:-1], h*1./len(bonds), color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L'); ylabel('pdf/#bonds')
plot([3,30], np.array([3, 30])**(1./3), 'k-')
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    h,bins = np.histogram(lengths, bins = np.arange(-1, max(lengths)+1))
    plot(bins[:-1], h, color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L'); ylabel('pdf')
plot([3,30], np.array([3, 30])**(-1./3), 'k-')
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    h,bins = np.histogram(lengths, bins = np.arange(-1, max(lengths)+1))
    plot(bins[:-1], h, color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L'); ylabel('pdf')
plot([3,30], 1e3*np.array([3, 30])**(-1./3), 'k-')
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    h,bins = np.histogram(lengths, bins = np.arange(-1, max(lengths)+1))
    plot(bins[:-1], h, color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L'); ylabel('pdf')
plot([3,30], 1e3*np.array([3, 30])**(-1./3), 'k-')
plot([3,30], 1e3*np.array([3, 30])**(-3), 'k-')
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    h,bins = np.histogram(lengths, bins = np.arange(-1, max(lengths)+1))
    plot(bins[:-1], h, color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L'); ylabel('pdf')
plot([3,30], 1e3*np.array([3, 30])**(-1./3), 'k-')
plot([3,30], 1e3*np.array([3, 30])**(-3.), 'k-')
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    h,bins = np.histogram(lengths, bins = np.arange(-1, max(lengths)+1))
    plot(bins[:-1], h, color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L'); ylabel('pdf')
plot([3,30], 1e3*np.array([3, 30])**(-1./3), 'k-')
plot([3,30], 1e3*np.array([3, 30])**(-5.), 'k-')
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    h,bins = np.histogram(lengths, bins = np.arange(-1, max(lengths)+1))
    plot(bins[:-1], h, color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L'); ylabel('pdf')
plot([3,30], 1e3*np.array([3, 30])**(-1./3), 'k-')
plot([3,30], 1e5*np.array([3, 30])**(-5.), 'k-')
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    h,bins = np.histogram(lengths, bins = np.arange(-1, max(lengths)+1))
    plot(bins[:-1], h, color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L'); ylabel('pdf')
plot([2,20], 1e3*np.array([2, 20])**(-1./3), 'k-')
plot([2,20], 1e5*np.array([2, 20])**(-5.), 'k-')
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    h,bins = np.histogram(lengths, bins = np.arange(-1, max(lengths)+1))
    plot(bins[:-1], h, color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L'); ylabel('pdf')
plot([2,20], 1e2*np.array([2, 20])**(-1./3), 'k-')
plot([2,20], 1e5*np.array([2, 20])**(-5.), 'k-')
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    h,bins = np.histogram(lengths, bins = np.arange(-1, max(lengths)+1))
    plot(bins[:-1], h, color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L'); ylabel('pdf')
plot([2,20], 1e2*np.array([2, 20])**(-1./2), 'k-')
plot([2,20], 1e5*np.array([2, 20])**(-5.), 'k-')
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    h,bins = np.histogram(lengths, bins = np.arange(-1, max(lengths)+1))
    plot(bins[:-1], h, color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L'); ylabel('pdf')
plot([2,20], 1e2*np.array([2, 20])**(-3./2), 'k-')
plot([2,20], 1e5*np.array([2, 20])**(-5.), 'k-')
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    h,bins = np.histogram(lengths, bins = np.arange(-1, max(lengths)+1))
    plot(bins[:-1], h, color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L'); ylabel('pdf')
plot([2,20], 5e2*np.array([2, 20])**(-3./2), 'k-')
plot([2,20], 1e5*np.array([2, 20])**(-5.), 'k-')
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    h,bins = np.histogram(lengths, bins = np.arange(-1, max(lengths)+1))
    plot(bins[:-1], h, color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L'); ylabel('pdf')
plot([2,20], 5e2*np.array([2, 20])**(-3./2), 'k-')
plot([2,20], 1e4*np.array([2, 20])**(-5.), 'k-')
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    h,bins = np.histogram(lengths, bins = np.arange(-1, max(lengths)+1))
    plot(bins[:-1], h, color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L'); ylabel('pdf')
plot([2,20], 5e2*np.array([2, 20])**(-3./2), 'k-')
plot([2,10], 1e4*np.array([2, 10])**(-5.), 'k-')
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    h,bins = np.histogram(lengths, bins = np.arange(-1, max(lengths)+1))
    plot(bins[:-1], h, color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L'); ylabel('pdf')
plot([2,20], 5e2*np.array([2, 20])**(-1.), 'k-')
plot([2,10], 1e4*np.array([2, 10])**(-5.), 'k-')
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    h,bins = np.histogram(lengths, bins = np.arange(-1, max(lengths)+1))
    plot(bins[:-1], h, color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L'); ylabel('pdf')
plot([2,20], 5e2*np.array([2, 20])**(-3./2), 'k-')
plot([2,10], 1e4*np.array([2, 10])**(-5.), 'k-')
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    h,bins = np.histogram(lengths, bins = np.arange(-1, max(lengths)+1))
    plot(bins[:-1]-1, h, color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([2,20], 5e2*np.array([2, 20])**(-3./2), 'k-')
plot([2,10], 1e4*np.array([2, 10])**(-5.), 'k-')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    h,bins = np.histogram(lengths, bins = np.arange(-1, max(lengths)+1))
    plot(bins[:-1]-1, h, color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,20], 5e2*np.array([1, 20])**(-3./2), 'k-')
plot([1,10], 1e3*np.array([1, 10])**(-5.), 'k-')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    h,bins = np.histogram(lengths, bins = np.arange(-1, max(lengths)+1))
    plot(bins[:-1]-1, h, color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,20], 5e2*np.array([1, 20])**(-3./2), 'k-')
plot([1,10], 5e2*np.array([1, 10])**(-3.), 'k-')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    h,bins = np.histogram(lengths, bins = np.arange(-1, max(lengths)+1))
    plot(bins[:-1]-1, h, color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,20], 5e2*np.array([1, 20])**(-3./2), 'k-')
plot([1,10], 5e2*np.array([1, 10])**(-4.), 'k-')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    h,bins = np.histogram(lengths, bins = np.arange(-1, max(lengths)+1))
    plot(bins[:-1]-1, h, color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,20], 5e2*np.array([1, 20])**(-3./2), 'k-')
plot([1,4], 5e2*np.array([1, 4])**(-4.), 'k-')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    h,bins = np.histogram(lengths, bins = np.arange(-1, max(lengths)+1))
    plot(bins[:-1]-1, h, color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,20], 5e2*np.array([1, 20])**(-3./2), 'k-')
plot([1,5], 5e2*np.array([1, 5])**(-4.), 'k-')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    h,bins = np.histogram(lengths, bins = np.arange(-1, max(lengths)+1))
    plot(bins[:-1]-1, h, color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,20], 5e2*np.array([1, 20])**(2.), 'k-')
plot([1,5], 5e2*np.array([1, 5])**(-4.), 'k-')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    h,bins = np.histogram(lengths, bins = np.arange(-1, max(lengths)+1))
    plot(bins[:-1]-1, h, color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,20], 5e2*np.array([1, 20])**(-2.), 'k-')
plot([1,5], 5e2*np.array([1, 5])**(-4.), 'k-')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    h,bins = np.histogram(lengths, bins = np.arange(-1, max(lengths)+1))
    plot(bins[:-1]-1, h, color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,20], 5e2*np.array([1, 20])**(-1.5), 'k-')
plot([1,5], 5e2*np.array([1, 5])**(-4.), 'k-')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    h,bins = np.histogram(lengths, bins = np.arange(-1, max(lengths)+1))
    plot(bins[:-1]-1, h*(bins[:-1]-1)**4, color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,20], 5e2*np.array([1, 20])**(-1.5), 'k-')
plot([1,5], 5e2*np.array([1, 5])**(-4.), 'k-')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    h,bins = np.histogram(lengths, bins = np.arange(-1, max(lengths)+1))
    plot(bins[:-1]-1, h*(bins[:-1]-1)**3, color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,20], 5e2*np.array([1, 20])**(-1.5), 'k-')
plot([1,5], 5e2*np.array([1, 5])**(-4.), 'k-')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    h,bins = np.histogram(lengths, bins = np.arange(-1, max(lengths)+1))
    plot(bins[:-1]-1, h*(bins[:-1]-1)**4, color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,20], 5e2*np.array([1, 20])**(-1.5), 'k-')
plot([1,5], 5e2*np.array([1, 5])**(-4.), 'k-')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    h,bins = np.histogram(lengths, bins = np.arange(-1, max(lengths)+1))
    plot(bins[:-1]-1, h*(bins[:-1]-1)**4, color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
#plot([1,20], 5e2*np.array([1, 20])**(-1.5), 'k-')
#plot([1,5], 5e2*np.array([1, 5])**(-4.), 'k-')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    h,bins = np.histogram(lengths, bins = np.arange(-1, max(lengths)+1))
    plot(bins[:-1]-1, h*(bins[:-1]-1)**4, color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,20], 5e2*np.array([1, 20])**(1.5), 'k-')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    h,bins = np.histogram(lengths, bins = np.arange(-1, max(lengths)+1))
    plot(bins[:-1]-1, h*(bins[:-1]-1)**4, color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,20], np.array([1, 20])**(2), 'k-')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    h,bins = np.histogram(lengths, bins = np.arange(-1, max(lengths)+1))
    plot(bins[:-1]-1, h*(bins[:-1]-1)**4, color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,20], 10*np.array([1, 20])**(2), 'k-')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    h,bins = np.histogram(lengths, bins = np.arange(-1, max(lengths)+1))
    plot(bins[:-1]-1, h*(bins[:-1]-1)**4, color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,20], 10*np.array([1, 20])**(3), 'k-')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    ml = max(lengths)+1
    ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
    h,bins = np.histogram(lengths, bins = ls)
    plot(bins[:-1]-1, h, color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,20], 5e2*np.array([1, 20])**(-1.5), 'k-')
plot([1,5], 5e2*np.array([1, 5])**(-4.), 'k-')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    ml = max(lengths)+1
    ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
    h,bins = np.histogram(lengths, bins = ls)
    h = h*1.*np.diff(bins)
    plot(bins[:-1]-1, h, color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,20], 5e2*np.array([1, 20])**(-1.5), 'k-')
plot([1,5], 5e2*np.array([1, 5])**(-4.), 'k-')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    ml = maxlength+1
    ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
    h,bins = np.histogram(lengths, bins = ls)
    h = h*1.*np.diff(bins)
    plot(bins[:-1]-1, h, color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,20], 5e2*np.array([1, 20])**(-1.5), 'k-')
plot([1,5], 5e2*np.array([1, 5])**(-4.), 'k-')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    ml = maxlength+1
    ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
    h,bins = np.histogram(lengths, bins = ls)
    h = h*1./np.diff(bins)
    plot(bins[:-1]-1, h, color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,20], 5e2*np.array([1, 20])**(-1.5), 'k-')
plot([1,5], 5e2*np.array([1, 5])**(-4.), 'k-')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    ml = maxlength+1
    ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
    h,bins = np.histogram(lengths, bins = ls)
    h = h*1./np.diff(bins)
    plot(bins[:-1]-1, h, color=cm.gist_earth(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,20], 5e2*np.array([1, 20])**(-1.5), 'k-')
plot([1,5], 5e2*np.array([1, 5])**(-4.), 'k-')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    ml = maxlength+1
    ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
    h,bins = np.histogram(lengths, bins = ls)
    h = h*1./np.diff(bins)
    plot(bins[:-1]-1, h, color=cm.gist_earth(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,20], 5e2*np.array([1, 20])**(-1.5), 'k:')
plot([1,5], 5e2*np.array([1, 5])**(-4.), 'k:')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    ml = maxlength+1
    ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
    h,bins = np.histogram(lengths, bins = ls)
    h = h*1./np.diff(bins)
    plot(bins[:-1]-1, h, color=cm.gist_earth(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,20], 5e2*np.array([1, 20])**(-1.5), 'k:')
plot([1,5], 4e2*np.array([1, 5])**(-3.), 'k:')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    ml = maxlength+1
    ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
    h,bins = np.histogram(lengths, bins = ls)
    h = h*1./np.diff(bins)
    plot(bins[:-1]-1, h, color=cm.gist_earth(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,20], 5e2*np.array([1, 20])**(-1.5), 'k:')
plot([1,5], 3e2*np.array([1, 5])**(-3.), 'k:')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    ml = maxlength+1
    ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
    h,bins = np.histogram(lengths, bins = ls)
    h = h*1./np.diff(bins)
    plot(bins[:-1]-1, h, color=cm.gist_earth(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,20], 5e2*np.array([1, 20])**(-1.5), 'k:')
plot([1,5], 3e2*np.array([1, 5])**(-3.5), 'k:')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    ml = maxlength+1
    ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
    h,bins = np.histogram(lengths, bins = ls)
    h = h*1./np.diff(bins)
    plot(bins[:-1]-1, h, color=cm.gist_earth(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,20], 5e2*np.array([1, 20])**(-1.5), 'k:')
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    ml = maxlength+1
    ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
    h,bins = np.histogram(lengths, bins = ls)
    h = h*1./np.diff(bins)
    plot(bins[:-1]-1, h, color=cm.gist_earth(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 5e2*np.array([1, 40])**(-1.5), 'k:')
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(bins)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 5e2*np.array([1, 40])**(-1.5), 'k:')
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(bins)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 5e2*np.array([1, 40])**(-1.5), 'k:')
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
for i, t in enumerate(ts):
    lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%t, int)
    h,bins = np.histogram(lengths, bins = np.arange(-1, max(lengths)+1))
    plot(bins[:-1]-1, h*(bins[:-1]-1)**3.5, color=cm.autumn(i*1./len(ts)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,20], 10*np.array([1, 20])**(3), 'k-')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(bins)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 1e2*np.array([1, 40])**(-0.5), 'k:')
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(bins)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 1e2*np.array([1, 40])**(-0.5), 'k:')
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(ls)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 1e2*np.array([1, 40])**(-0.5), 'k:')
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(ls)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 1e2*np.array([1, 40])**(-1.), 'k:')
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(ls)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 2e2*np.array([1, 40])**(-1.), 'k:')
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(ls)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 1.5e2*np.array([1, 40])**(-1.), 'k:')
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(ls)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 2e2*np.array([1, 40])**(-1.), 'k:')
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(ls)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 2e2*np.array([1, 40])**(-1.), 'k:')
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
plot([1,40], 2e1*np.array([1, 40])**(-1.), 'k:')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(ls)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 2e2*np.array([1, 40])**(-1.), 'k:')
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
plot([1,40], 2e1*np.array([1, 40])**(-2.), 'k:')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(ls)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 2e2*np.array([1, 40])**(-1.), 'k:')
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
plot([1,40], 2e1*np.array([1, 40])**(-1.5), 'k:')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(ls)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 2e2*np.array([1, 40])**(-1.), 'k:')
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
plot([1,40], 5e1*np.array([1, 40])**(-1.5), 'k:')
L=5
t0 = 68
for na in ['', 'bulk']:    
    lengths = [
        np.loadtxt(name, int) 
        for t, name in x.enum('_broken'+na, 'length') 
        if t>0]
    maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
    his = np.array([
        np.histogram(l, np.arange(-1, maxlength+1))[0]
        for l in lengths])
    r = his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(li[:-1], [r[i:j].mean() for i,j in zip(li, li[1:])], 'o-', label=lab)[0]
xscale('log')
yscale('log')
#ylim(1e-6,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(ls)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 2e2*np.array([1, 40])**(-1.), 'k:')
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
plot([1,40], 5e1*np.array([1, 40])**(-1.5), 'k:')
plot([10,70], 1e5*np.array([10, 70])**(-4), 'k:')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(ls)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 2e2*np.array([1, 40])**(-1.), 'k:')
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
plot([1,40], 5e1*np.array([1, 40])**(-1.5), 'k:')
plot([10,70], 1e5*np.array([10, 70])**(-4.), 'k:')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(ls)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 2e2*np.array([1, 40])**(-1.), 'k:')
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
plot([1,40], 5e1*np.array([1, 40])**(-1.5), 'k:')
plot([10,70], 1e6*np.array([10, 70])**(-4.), 'k:')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(ls)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 2e2*np.array([1, 40])**(-1.), 'k:')
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
plot([1,40], 5e1*np.array([1, 40])**(-1.5), 'k:')
plot([10,70], 1e6*np.array([10, 70])**(-5.), 'k:')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(ls)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 2e2*np.array([1, 40])**(-1.), 'k:')
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
plot([1,40], 5e1*np.array([1, 40])**(-1.5), 'k:')
plot([10,70], 1e7*np.array([10, 70])**(-5.), 'k:')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(ls)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 2e2*np.array([1, 40])**(-1.), 'k:')
text(3, 150, "-1")
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
plot([1,40], 5e1*np.array([1, 40])**(-1.5), 'k:')
plot([10,70], 1e7*np.array([10, 70])**(-5.), 'k:')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(ls)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 2e2*np.array([1, 40])**(-1.), 'k:')
text(3, 100, "-1")
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
plot([1,40], 5e1*np.array([1, 40])**(-1.5), 'k:')
plot([10,70], 1e7*np.array([10, 70])**(-5.), 'k:')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(ls)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 2e2*np.array([1, 40])**(-1.), 'k:')
text(3, 100, "-1")
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
text(1.5, 20, "-3.5")
plot([1,40], 5e1*np.array([1, 40])**(-1.5), 'k:')
plot([10,70], 1e7*np.array([10, 70])**(-5.), 'k:')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(ls)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 2e2*np.array([1, 40])**(-1.), 'k:')
text(3, 100, "-1")
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
text(1.2, 300, "-3.5")
plot([1,40], 5e1*np.array([1, 40])**(-1.5), 'k:')
text(1.2, 20, "-1.5")
plot([10,70], 1e7*np.array([10, 70])**(-5.), 'k:')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(ls)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 2e2*np.array([1, 40])**(-1.), 'k:')
text(3, 100, "-1")
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
text(5, 1, "-3.5")
plot([1,40], 5e1*np.array([1, 40])**(-1.5), 'k:')
text(1.2, 20, "-1.5")
plot([10,70], 1e7*np.array([10, 70])**(-5.), 'k:')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(ls)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 2e2*np.array([1, 40])**(-1.), 'k:')
text(3, 100, "-1")
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
text(4, 1, "-3.5")
plot([1,40], 5e1*np.array([1, 40])**(-1.5), 'k:')
text(1.2, 20, "-1.5")
plot([10,70], 1e7*np.array([10, 70])**(-5.), 'k:')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(ls)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 2e2*np.array([1, 40])**(-1.), 'k:')
text(3, 100, "-1")
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
text(4, 1, "-3.5")
plot([1,40], 5e1*np.array([1, 40])**(-1.5), 'k:')
text(1.2, 20, "-1.5")
plot([10,70], 1e7*np.array([10, 70])**(-5.), 'k:')
text(60, 5e-3, "-5")
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(ls)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 2e2*np.array([1, 40])**(-1.), 'k:')
text(3, 100, "-1")
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
text(4, 1, "-3.5")
plot([1,40], 5e1*np.array([1, 40])**(-1.5), 'k:')
text(1.2, 20, "-1.5")
plot([10,70], 1e7*np.array([10, 70])**(-5.), 'k:')
text(60, 5e-3, "-5")
plot([10,70], 1e7*np.array([10, 70])**(-6.), 'k:')
text(10, 5e-3, "-6")
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(ls)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 2e2*np.array([1, 40])**(-1.), 'k:')
text(3, 100, "-1")
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
text(4, 1, "-3.5")
plot([1,40], 5e1*np.array([1, 40])**(-1.5), 'k:')
text(1.2, 20, "-1.5")
plot([10,70], 1e7*np.array([10, 70])**(-5.), 'k:')
text(60, 5e-3, "-5")
plot([5,10], 1e5*np.array([5, 10])**(-6.), 'k:')
text(10, 5e-3, "-6")
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(ls)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 2e2*np.array([1, 40])**(-1.), 'k:')
text(3, 100, "-1")
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
text(4, 1, "-3.5")
plot([1,40], 5e1*np.array([1, 40])**(-1.5), 'k:')
text(1.2, 20, "-1.5")
plot([10,70], 1e7*np.array([10, 70])**(-5.), 'k:')
text(60, 5e-3, "-5")
plot([5,10], 1e5*np.array([5, 10])**(-7.), 'k:')
text(10, 5e-3, "-6")
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(ls)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 2e2*np.array([1, 40])**(-1.), 'k:')
text(3, 100, "-1")
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
text(4, 1, "-3.5")
plot([1,40], 5e1*np.array([1, 40])**(-1.5), 'k:')
text(1.2, 20, "-1.5")
plot([10,70], 1e7*np.array([10, 70])**(-5.), 'k:')
text(60, 5e-3, "-5")
plot([5,10], 1e6*np.array([5, 10])**(-8.), 'k:')
text(10, 5e-3, "-6")
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(ls)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 2e2*np.array([1, 40])**(-1.), 'k:')
text(3, 100, "-1")
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
text(4, 1, "-3.5")
plot([1,40], 5e1*np.array([1, 40])**(-1.5), 'k:')
text(1.2, 20, "-1.5")
plot([10,70], 1e7*np.array([10, 70])**(-5.), 'k:')
text(60, 5e-3, "-5")
plot([5,10], 1e6*np.array([5, 10])**(-9.), 'k:')
text(10, 5e-3, "-6")
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(ls)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 2e2*np.array([1, 40])**(-1.), 'k:')
text(3, 100, "-1")
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
text(4, 1, "-3.5")
plot([1,40], 5e1*np.array([1, 40])**(-1.5), 'k:')
text(1.2, 20, "-1.5")
plot([10,70], 1e7*np.array([10, 70])**(-5.), 'k:')
text(60, 5e-3, "-5")
plot([5,10], 1e8*np.array([5, 10])**(-9.), 'k:')
text(10, 5e-3, "-9")
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(ls)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 2e2*np.array([1, 40])**(-1.), 'k:')
text(3, 100, "-1")
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
text(4, 1, "-3.5")
plot([1,40], 5e1*np.array([1, 40])**(-1.5), 'k:')
text(1.2, 20, "-1.5")
plot([10,70], 1e7*np.array([10, 70])**(-5.), 'k:')
text(60, 5e-3, "-5")
plot([8,15], 1e8*np.array([8, 15])**(-9.), 'k:')
text(10, 5e-3, "-9")
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_brokenbulk', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(ls)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 2e2*np.array([1, 40])**(-1.), 'k:')
text(3, 100, "-1")
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
text(4, 1, "-3.5")
plot([1,40], 5e1*np.array([1, 40])**(-1.5), 'k:')
text(1.2, 20, "-1.5")
plot([10,70], 1e7*np.array([10, 70])**(-5.), 'k:')
text(60, 5e-3, "-5")
plot([8,15], 1e8*np.array([8, 15])**(-9.), 'k:')
text(10, 5e-3, "-9")
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/breaking_L_t.pdf')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_broken', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(ls)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 2e2*np.array([1, 40])**(-1.), 'k:')
text(3, 100, "-1")
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
text(4, 1, "-3.5")
plot([1,40], 5e1*np.array([1, 40])**(-1.5), 'k:')
text(1.2, 20, "-1.5")
plot([10,70], 1e7*np.array([10, 70])**(-5.), 'k:')
text(60, 5e-3, "-5")
plot([8,15], 1e8*np.array([8, 15])**(-9.), 'k:')
text(10, 5e-3, "-9")
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/breaking_L_t.pdf')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_broken', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(ls)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 2e2*np.array([1, 40])**(-1.), 'k:')
text(3, 100, "-1")
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
text(4, 1, "-3.5")
plot([1,40], 5e1*np.array([1, 40])**(-1.5), 'k:')
text(1.2, 20, "-1.5")
plot([10,70], 1e5*np.array([10, 70])**(-3.), 'k:')
text(60, 5e-3, "-5")
plot([8,15], 1e8*np.array([8, 15])**(-9.), 'k:')
text(10, 5e-3, "-9")
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/breaking_L_t.pdf')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_broken', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(ls)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 2e2*np.array([1, 40])**(-1.), 'k:')
text(3, 100, "-1")
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
text(4, 1, "-3.5")
plot([1,40], 5e1*np.array([1, 40])**(-1.5), 'k:')
text(1.2, 20, "-1.5")
plot([10,70], 1e4*np.array([10, 70])**(-3.), 'k:')
text(60, 5e-3, "-5")
plot([8,15], 1e8*np.array([8, 15])**(-9.), 'k:')
text(10, 5e-3, "-9")
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/breaking_L_t.pdf')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_broken', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(ls)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 2e2*np.array([1, 40])**(-1.), 'k:')
text(3, 100, "-1")
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
text(4, 1, "-3.5")
plot([1,40], 5e1*np.array([1, 40])**(-1.5), 'k:')
text(1.2, 20, "-1.5")
plot([10,70], 5e4*np.array([10, 70])**(-3.), 'k:')
text(60, 5e-3, "-5")
plot([8,15], 1e8*np.array([8, 15])**(-9.), 'k:')
text(10, 5e-3, "-9")
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/breaking_L_t.pdf')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_broken', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(ls)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 2e2*np.array([1, 40])**(-1.), 'k:')
text(3, 100, "-1")
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
text(4, 1, "-3.5")
plot([1,40], 5e1*np.array([1, 40])**(-1.5), 'k:')
text(1.2, 20, "-1.5")
plot([10,70], 5e4*np.array([10, 70])**(-3.), 'k:')
text(60, 5e-1, "-3")
plot([8,15], 1e5*np.array([8, 15])**(-6.), 'k:')
text(10, 5e-3, "-9")
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/breaking_L_t.pdf')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_broken', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(ls)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 2e2*np.array([1, 40])**(-1.), 'k:')
text(3, 100, "-1")
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
text(4, 1, "-3.5")
plot([1,40], 5e1*np.array([1, 40])**(-1.5), 'k:')
text(1.2, 20, "-1.5")
plot([10,70], 5e4*np.array([10, 70])**(-3.), 'k:')
text(60, 5e-1, "-3")
plot([8,15], 1e5*np.array([8, 15])**(-6.), 'k:')
text(10, 5e-2, "-6")
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/breaking_L_t.pdf')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_broken', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(ls)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 2e2*np.array([1, 40])**(-1.), 'k:')
text(3, 100, "-1")
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
text(4, 1, "-3.5")
plot([1,40], 5e1*np.array([1, 40])**(-1.5), 'k:')
text(1.2, 20, "-1.5")
plot([10,70], 5e4*np.array([10, 70])**(-3.), 'k:')
text(60, 5e-1, "-3")
plot([8,15], 1e5*np.array([8, 15])**(-6.), 'k:')
text(10, 3e-2, "-6")
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/breaking_L_t.pdf')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_broken', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(ls)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 2e2*np.array([1, 40])**(-1.), 'k:')
text(3, 100, "-1")
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
text(4, 1, "-3.5")
plot([1,40], 5e1*np.array([1, 40])**(-1), 'k:')
text(1.2, 20, "-1.5")
plot([10,70], 5e4*np.array([10, 70])**(-3.), 'k:')
text(60, 5e-1, "-3")
plot([8,15], 1e5*np.array([8, 15])**(-6.), 'k:')
text(10, 3e-2, "-6")
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/breaking_L_t.pdf')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_broken', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(ls)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 2e2*np.array([1, 40])**(-1.), 'k:')
text(3, 100, "-1")
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
text(4, 1, "-3.5")
plot([1,40], 5e1*np.array([1, 40])**(-1.), 'k:')
text(1.2, 20, "-1.5")
plot([10,70], 5e4*np.array([10, 70])**(-3.), 'k:')
text(60, 5e-1, "-3")
plot([8,15], 1e5*np.array([8, 15])**(-6.), 'k:')
text(10, 3e-2, "-6")
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/breaking_L_t.pdf')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_broken', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
    h = h*1./np.diff(ls)/(ts[i+1]-t)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
plot([1,40], 2e2*np.array([1, 40])**(-1.), 'k:')
text(3, 100, "-1")
plot([1,5], 4e2*np.array([1, 5])**(-3.5), 'k:')
text(4, 1, "-3.5")
plot([1,40], 5e1*np.array([1, 40])**(-1.5), 'k:')
text(1.2, 20, "-1.5")
plot([10,70], 5e4*np.array([10, 70])**(-3.), 'k:')
text(60, 5e-1, "-3")
plot([8,15], 1e5*np.array([8, 15])**(-6.), 'k:')
text(10, 3e-2, "-6")
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/breaking_L_t.pdf')
L=5
t0 = 68
for na in ['', 'bulk']:    
    lengths = [
        np.loadtxt(name, int) 
        for t, name in x.enum('_broken'+na, 'length') 
        if t>0]
    maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
    his = np.array([
        np.histogram(l, np.arange(-1, maxlength+1))[0]
        for l in lengths])
    r = his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(li[:-1], [r[i:j].mean() for i,j in zip(li, li[1:])], 'o-', label=lab)[0]
xscale('log')
yscale('log')
#ylim(1e-6,1e-1)
#xlim(1,1e3)
xlabel(r'$t/\tau_B$')
ylabel(r'strand rupture per particle')
print li
x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
lengths = [
    np.loadtxt(name, int) 
    for t, name in x.enum('_broken', 'length') 
    if t>0]
maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
his = np.array([
    np.histogram(l, np.arange(-1, maxlength+1))[0]
    for l in lengths])
plot(np.arange(maxlength-1), his[t0:t0+9,2:].mean(0))
plot(np.arange(maxlength-1), his[t0+9:t0+40,2:].mean(0))
plot(np.arange(maxlength-1), his[t0+40:,2:].mean(0))
xscale('log');yscale('log')
xlabel('L')
ylabel('pdf')
figure(figsize(8,6))
x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
lengths = [
    np.loadtxt(name, int) 
    for t, name in x.enum('_broken', 'length') 
    if t>0]
maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
his = np.array([
    np.histogram(l, np.arange(-1, maxlength+1))[0]
    for l in lengths])
plot(np.arange(maxlength-1), his[t0:t0+9,2:].mean(0))
plot(np.arange(maxlength-1), his[t0+9:t0+40,2:].mean(0))
plot(np.arange(maxlength-1), his[t0+40:,2:].mean(0))
xscale('log');yscale('log')
xlabel('L')
ylabel('pdf')
figure(figsize(8,6))
x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
lengths = [
    np.loadtxt(name, int) 
    for t, name in x.enum('_broken', 'length') 
    if t>0]
maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
his = np.array([
    np.histogram(l, np.arange(-1, maxlength+1))[0]
    for l in lengths])
plot(np.arange(maxlength-1), his[t0:t0+9,2:].mean(0))
plot(np.arange(maxlength-1), his[t0+9:t0+40,2:].mean(0))
plot(np.arange(maxlength-1), his[t0+40:,2:].mean(0))
xscale('log');yscale('log')
xlabel('L')
ylabel('pdf')
xlim(1,100)
lengths = np.loadtxt(x.get_format_string('_broken', 'length')%111, int)
h = np.histogram(lengths, bins = ls)[0]
plot(ls[:-1]-1, h)
lengths = np.loadtxt(x.get_format_string('_broken', 'length')%111, int)
h = np.histogram(lengths, bins = ls)[0]
plot(ls[:-1]-1, h)
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
lengths = np.loadtxt(x.get_format_string('_broken', 'length')%(t0+9), int)
h = np.histogram(lengths, bins = ls)[0]
plot(ls[:-1]-1, h)
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
for t in [0, 10, 136]:
    lengths = np.loadtxt(x.get_format_string('_broken', 'length')%(t0+t), int)
    h = np.histogram(lengths, bins = ls)[0]
    plot(ls[:-1]-1, h)
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
for t in [0, 10, 40, 136]:
    lengths = np.loadtxt(x.get_format_string('_broken', 'length')%(t0+t), int)
    h = np.histogram(lengths, bins = ls)[0]
    plot(ls[:-1]-1, h)
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
for t in [0, 10, 40]:
    lengths = np.loadtxt(x.get_format_string('_broken', 'length')%(t0+t), int)
    h = np.histogram(lengths, bins = ls)[0]
    plot(ls[:-1]-1, h)
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
for t in [0, 10, 40]:
    lengths = np.loadtxt(x.get_format_string('_broken', 'length')%(t0+t), int)
    bonds = np.loadtxt(x.get_format_string(exxt='bonds'), dtype=int)%(t0+t), int)
    h = np.histogram(lengths, bins = ls)[0]
    plot(ls[:-1]-1, h*1./len(bonds))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
for t in [0, 10, 40]:
    lengths = np.loadtxt(x.get_format_string('_broken', 'length')%(t0+t), int)
    bonds = np.loadtxt(x.get_format_string(exxt='bonds')%(t0+t), int)
    h = np.histogram(lengths, bins = ls)[0]
    plot(ls[:-1]-1, h*1./len(bonds))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
for t in [0, 10, 40]:
    lengths = np.loadtxt(x.get_format_string('_broken', 'length')%(t0+t), int)
    bonds = np.loadtxt(x.get_format_string(ext='bonds')%(t0+t), int)
    h = np.histogram(lengths, bins = ls)[0]
    plot(ls[:-1]-1, h*1./len(bonds))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    nbbonds = 0
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_broken', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
        nbbonds += sum([1 for line in open(x.get_format_string(ext='bonds')%tt)])
    h = h*1./np.diff(ls)/(nbbonds)
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
figure(figsize=(8,6))
ts = t0 +1 + np.unique(np.logspace(0, np.log10(x.size-t0-1), 10*np.log10(x.size-t0-1)).astype(int))
ml = maxlength+1
ls = np.unique(np.logspace(0, np.log10(ml), 10*np.log10(ml)).astype(int))
for i, t in enumerate(ts[:-1]):
    h = np.zeros(len(ls)-1, int)
    nbbonds = 0
    for tt in range(t, ts[i+1]):
        lengths = np.loadtxt(x.get_format_string('_broken', 'length')%tt, int)
        h += np.histogram(lengths, bins = ls)[0]
        nbbonds += sum([1 for line in open(x.get_format_string(ext='bonds')%tt)])
    h = h*1./np.diff(ls)/nbbonds
    plot(ls[:-1]-1, h, color=cm.gist_earth(i*1./(len(ts)-1)))
xscale('log'); yscale('log')
xlabel('L-1'); ylabel('pdf')
x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/2_ageing/1415_ageing.traj')
bonds0 = np.atleast_2d(np.loadtxt(x.get_format_string(ext='bonds')%0, dtype=int))
bonds1 = np.atleast_2d(np.loadtxt(x.get_format_string(ext='bonds')%1, dtype=int))
p2tr0 = np.loadtxt(x.get_format_string(ext='p2tr')%0, dtype=int)
p2tr1 = np.loadtxt(x.get_format_string(ext='p2tr')%1, dtype=int)
trbonds = set([(a,b) for a,b in np.sort(p2tr0[bonds0], 1)]) - set([(a,b) for a,b in np.sort(p2tr1[bonds1], axis=1)])
#graph of the bonds between trajectories at t+dt
g = nx.Graph()
g.add_nodes_from(p2tr1)
g.add_edges_from(p2tr1[bonds1])
trajs = np.unique([[trpath[0], trpath[-1]] for trpath in shortest_path(g, trbonds)])
x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/2_ageing/1415_ageing.traj')
bonds0 = np.atleast_2d(np.loadtxt(x.get_format_string(ext='bonds')%0, dtype=int))
bonds1 = np.atleast_2d(np.loadtxt(x.get_format_string(ext='bonds')%1, dtype=int))
p2tr0 = np.loadtxt(x.get_format_string(ext='p2tr')%0, dtype=int)
p2tr1 = np.loadtxt(x.get_format_string(ext='p2tr')%1, dtype=int)
trbonds = set([(a,b) for a,b in np.sort(p2tr0[bonds0], 1)]) - set([(a,b) for a,b in np.sort(p2tr1[bonds1], axis=1)])
#graph of the bonds between trajectories at t+dt
g = nx.Graph()
g.add_nodes_from(p2tr1)
g.add_edges_from(p2tr1[bonds1])
trajs = np.array([[trpath[0], trpath[-1]] for trpath in shortest_path(g, trbonds)])
positions = [np.zeros((len(x.trajs[tr]), 3)) for tr in trajs] #do not take care of x.trajstart since we start at t=0
for t, name in x.enum():
    pos = np.loadtxt(name, skiprows=2)
    for itr, tr in enumerate(trajs):
        if t >= len(tr): continue
        positions[itr][t] = pos[x.trajs[tr][t]]
get_ipython().magic(u'debug ')
positions = [[np.zeros((len(x.trajs[tr]), 3)) for tr in trcouple] for trclouple in trajs] #do not take care of x.trajstart since we start at t=0
for t, name in x.enum():
    pos = np.loadtxt(name, skiprows=2)
    for itrcouple, trcouple in enumerate(trajs):
        for itr, tr in trcouple:
            if t >= len(tr): continue
            positions[itrcouple][itr][t] = pos[x.trajs[tr][t]]
positions = [[np.zeros((len(x.trajs[tr]), 3)) for tr in trcouple] for trcouple in trajs] #do not take care of x.trajstart since we start at t=0
for t, name in x.enum():
    pos = np.loadtxt(name, skiprows=2)
    for itrcouple, trcouple in enumerate(trajs):
        for itr, tr in trcouple:
            if t >= len(tr): continue
            positions[itrcouple][itr][t] = pos[x.trajs[tr][t]]
positions = [[np.zeros((len(x.trajs[tr]), 3)) for tr in trcouple] for trcouple in trajs] #do not take care of x.trajstart since we start at t=0
for t, name in x.enum():
    pos = np.loadtxt(name, skiprows=2)
    for itrcouple, trcouple in enumerate(trajs):
        for itr, tr in enumerate(trcouple):
            if t >= len(tr): continue
            positions[itrcouple][itr][t] = pos[x.trajs[tr][t]]
positions = [[np.zeros((len(x.trajs[tr]), 3)) for tr in trcouple] for trcouple in trajs] #do not take care of x.trajstart since we start at t=0
for t, name in x.enum():
    pos = np.loadtxt(name, skiprows=2)
    for itrcouple, trcouple in enumerate(trajs):
        for itr, tr in enumerate(trcouple):
            if t >= len(x.trajs[tr]): continue
            positions[itrcouple][itr][t] = pos[x.trajs[tr][t]]
print positions[0][0].shape
print positions[0][1].shape
np.argmax(map(len, shortest_path(g, trbonds)))
print positions[191][0].shape
print positions[191][1].shape
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(positions[191][0][:,0], positions[191][0][:,1], positions[191][0][:,2])
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(positions[191][0][:,0], positions[191][0][:,1], positions[191][0][:,2])
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for ps in positions[191]:
    ax.plot(ps[:,0], ps[0][:,1], ps[:,2])
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for ps in positions[191]:
    ax.plot(ps[:,0], ps[:,1], ps[:,2])
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for ps in positions[0]:
    ax.plot(ps[:,0], ps[:,1], ps[:,2])
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for ps in positions[0]:
    ax.plot(ps[:,0], ps[:,1], ps[:,2], color=cm.autumn(np.linspace(0,1,118)))
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for ps in positions[0]:
    ax.scatter(ps[:,0], ps[:,1], ps[:,2], color=cm.autumn(np.linspace(0,1,118)))
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for ps in positions[191]:
    ax.scatter(ps[:,0], ps[:,1], ps[:,2], color=cm.autumn(np.linspace(0,1,118)[:len(ps)]))
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for ps in positions[1]:
    ax.scatter(ps[:,0], ps[:,1], ps[:,2], color=cm.autumn(np.linspace(0,1,118)[:len(ps)]))
plot(np.sum((positions[0][0] - positions[0][1])**2, -1))
plot(np.sqrt(np.sum((positions[0][0] - positions[0][1])**2, -1)))
plot(np.sqrt(np.sum((positions[0][0] - positions[0][1])**2, -1))/10.)
plot(np.sqrt(np.sum((positions[191][0] - positions[191][1])**2, -1))/10.)
plot(np.sqrt(np.sum((positions[191][0][:44] - positions[191][1][:44])**2, -1))/10.)
len(positions)
for pp in positions:
    lm = max(map(len, pp))
    plot(np.sqrt(np.sum((pp[0][:lm] - pp[1][:lm])**2, -1))/10.)
get_ipython().magic(u'debug ')
for pp in positions:
    lm = min(map(len, pp))
    plot(np.sqrt(np.sum((pp[0][:lm] - pp[1][:lm])**2, -1))/10.)
tot = np.zeros((100, x.size))
nb = np.zeros((100, x.size), int)
for pp, path in zip(positions, shortest_path(g, trbonds)):
    L = len(path)-1
    lm = min(map(len, pp))
    tot[L,:lm] += np.sum((pp[0][:lm] - pp[1][:lm])**2, -1)
    nb[L, :lm] += 1
meandist = np.sqrt(tot/nb)/10
tot = np.zeros((100, x.size))
nb = np.zeros((100, x.size), int)
for pp, path in zip(positions, shortest_path(g, trbonds)):
    L = len(path)-1
    lm = min(map(len, pp))
    tot[L,:lm] += np.sum((pp[0][:lm] - pp[1][:lm])**2, -1)
    nb[L, :lm] += 1
meandist = np.sqrt(tot/np.maximum(1, nb))/10
plot(np.sqrt(tot.sum(0)/np.maximum(1, nb.sum(0)))/10
plot(np.sqrt(tot.sum(0)/np.maximum(1, nb.sum(0)))/10)
plot(np.sqrt(tot.sum(0)/np.maximum(1, nb.sum(0)))/10)
plot(np.sqrt(tot[2:].sum(0)/np.maximum(1, nb[2:].sum(0)))/10)
plot(np.sqrt(tot[3:].sum(0)/np.maximum(1, nb[3:].sum(0)))/10)
plot(np.sqrt(tot.sum(0)/np.maximum(1, nb.sum(0)))/10)
for L in [3,13,23]:
    plot(np.sqrt(tot[L:].sum(0)/np.maximum(1, nb[L:].sum(0)))/10)
plot(np.sqrt(tot.sum(0)/np.maximum(1, nb.sum(0)))/10)
for L in [3,13,18]:
    plot(np.sqrt(tot[L:].sum(0)/np.maximum(1, nb[L:].sum(0)))/10)
np.max(map(len, shortest_path(g, trbonds)))
print np.argmax(map(len, shortest_path(g, trbonds)))
print np.max(map(len, shortest_path(g, trbonds)))
tot = np.zeros((21, x.size))
nb = np.zeros((21, x.size), int)
for pp, path in zip(positions, shortest_path(g, trbonds)):
    L = len(path)-1
    lm = min(map(len, pp))
    tot[L, :lm] += np.sum((pp[0][:lm] - pp[1][:lm])**2, -1)
    nb[L, :lm] += 1
meandist = np.sqrt(tot/np.maximum(1, nb))/10
plot(np.sqrt(tot.sum(0)/np.maximum(1, nb.sum(0)))/10)
for L in [3,5,15]:
    plot(np.sqrt(tot[L:].sum(0)/np.maximum(1, nb[L:].sum(0)))/10)
plot(np.sqrt(tot.sum(0)/np.maximum(1, nb.sum(0)))/10)
for L in [3,5,13]:
    plot(np.sqrt(tot[L:].sum(0)/np.maximum(1, nb[L:].sum(0)))/10)
plot(np.sqrt(tot.sum(0)/np.maximum(1, nb.sum(0)))/10)
for L in [3,5,10]:
    plot(np.sqrt(tot[L:].sum(0)/np.maximum(1, nb[L:].sum(0)))/10)
plot(np.sqrt(tot.sum(0)/np.maximum(1, nb.sum(0)))/10)
for L in [3,5,10]:
    plot(np.sqrt(tot[L:].sum(0)/np.maximum(1, nb[L:].sum(0)))/10)
xscale('log')
tot = np.zeros((21, x.size))
nb = np.zeros((21, x.size), int)
for pp, path in zip(positions, shortest_path(g, trbonds)):
    L = len(path)-1
    lm = min(map(len, pp))
    if lm<118: continue
    tot[L, :lm] += np.sum((pp[0][:lm] - pp[1][:lm])**2, -1)
    nb[L, :lm] += 1
meandist = np.sqrt(tot/np.maximum(1, nb))/10
plot(np.sqrt(tot.sum(0)/np.maximum(1, nb.sum(0)))/10)
for L in [3,5,10]:
    plot(np.sqrt(tot[L:].sum(0)/np.maximum(1, nb[L:].sum(0)))/10)
xscale('log')
tot = np.zeros((21, x.size))
nb = np.zeros((21, x.size), int)
for pp, path in zip(positions, shortest_path(g, trbonds)):
    L = len(path)-1
    lm = min(map(len, pp))
    depsq = np.sum((pp[0][:lm] - pp[1][:lm])**2, -1)
    if depsq[2]<depsq[1]: continue
    tot[L, :lm] += depsq
    nb[L, :lm] += 1
meandist = np.sqrt(tot/np.maximum(1, nb))/10
tot = np.zeros((21, x.size))
nb = np.zeros((21, x.size), int)
for pp, path in zip(positions, shortest_path(g, trbonds)):
    L = len(path)-1
    lm = min(map(len, pp))
    depsq = np.sum((pp[0][:lm] - pp[1][:lm])**2, -1)
    if lm<3 or depsq[2]<depsq[1]: continue
    tot[L, :lm] += depsq
    nb[L, :lm] += 1
meandist = np.sqrt(tot/np.maximum(1, nb))/10
plot(np.sqrt(tot.sum(0)/np.maximum(1, nb.sum(0)))/10)
for L in [3,5,10]:
    plot(np.sqrt(tot[L:].sum(0)/np.maximum(1, nb[L:].sum(0)))/10)
xscale('log')
plot(np.sqrt(tot.sum(0)/np.maximum(1, nb.sum(0)))/10)
for L in [3,5,6]:
    plot(np.sqrt(tot[L:].sum(0)/np.maximum(1, nb[L:].sum(0)))/10)
xscale('log')
plot(np.sqrt(tot.sum(0)/np.maximum(1, nb.sum(0)))/10)
for L in [3,4,5]:
    plot(np.sqrt(tot[L:].sum(0)/np.maximum(1, nb[L:].sum(0)))/10)
xscale('log')
plot(np.sqrt(tot.sum(0)/np.maximum(1, nb.sum(0)))/10)
for L in [3,4,5]:
    mdist = tot[L:].sum(0)/np.maximum(1, nb[L:].sum(0))
    plot(np.sqrt(mdist/mdist[0]))
xscale('log')
plot(np.sqrt(tot.sum(0)/np.maximum(1, nb.sum(0)))/10)
for L in [3,4,5]:
    plot(np.sqrt(tot[L:].sum(0)/np.maximum(1, nb[L:].sum(0)))/10)
xscale('log')
maxdepsq = 0
argmaxdep = 0
for i, (pp, path) in enumerate(zip(positions, shortest_path(g, trbonds))):
    L = len(path)-1
    lm = min(map(len, pp))
    if lm<118: continue
    depsq = np.sum((pp[0][:lm] - pp[1][:lm])**2, -1)
    if depsq[-1]>maxdepsq:
        maxdepsq = depsq[-1]
        argmaxdep = i
print argmaxdep, maxdepsq
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for ps in positions[189]:
    ax.scatter(ps[:,0], ps[:,1], ps[:,2], color=cm.autumn(np.linspace(0,1,118)[:len(ps)]))
np.sqrt(maxdepsq)/10
def broken_bonds_path(bonds0, bonds1, p2tr0, p2tr1):
    """Paths on graph at t1 between particles involved in a bond at t0 and no more at t1.
    
    bonds0, bonds1 are respectively the bonds at t0 and 1 in terms of position
    p2tr0, p2tr1 are respectively the position to trajectory relationship at t0 and t1
    
    generate paths in terms of trajectory indices"""
    #bonds (between trajectories) existing at t but no more at t+dt
    # = broken bonds + lost trajectories
    trbonds = set([(a,b) for a,b in np.sort(p2tr0[bonds0], 1)]) - set([(a,b) for a,b in np.sort(p2tr1[bonds1], axis=1)])

    #graph of the bonds between trajectories at t+dt
    g = nx.Graph()
    g.add_nodes_from(p2tr1)
    g.add_edges_from(p2tr1[bonds1])
    return shortest_path(g, trbonds)
def broken_bonds_path_p(bonds0, bonds1, p2tr0, p2tr1):
    """Paths on graph at t1 between particles involved in a bond at t0 and no more at t1.
    
    bonds0, bonds1 are respectively the bonds at t0 and 1 in terms of position
    p2tr0, p2tr1 are respectively the position to trajectory relationship at t0 and t1
    
    generate paths in terms of position indices in t1"""
    tr2p1 = dict((tr,p) for p, tr in enumerate(p2tr1))
    
    for trpath in broken_bonds_path(bonds0, bonds1, p2tr0, p2tr1):
        yield [tr2p1[tr] for tr in trpath]
sorted((1,0))
tuple(sorted((1,0)))
#look for all the trajectories that are ever bounded together and unbounded at the next time step
x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
pro = ProgressBar(x.size)
everbounded = set()
for t, name in x.enum(ext='bonds'):
    try:
        bonds1 = np.atleast_2d(np.loadtxt(name, dtype=int))
    except UserWarning:
        bonds1 = np.zeros([0,2], int)
    p2tr1 = np.loadtxt(x.get_format_string(ext='p2tr')%t, dtype=int)
    if t>0:
        everbounded |= set(tuple(sorted((path[0], path[-1]))) for path in broken_bonds_path(bonds0, bonds1, p2tr0, p2tr1))
    bonds0 = bonds1
    p2tr0 = p2tr1
    pro.animate(t)
len(x.trajstart), len(everbounded), len(everbounded)/float(len(x.trajstart))
len(x.nb_trajs), len(everbounded), len(everbounded)/float(len(x.trajstart))
x.nb_trajs, len(everbounded), len(everbounded)/float(x.nb_trajs))
x.nb_trajs, len(everbounded), len(everbounded)/float(x.nb_trajs)
everbounded_tr = np.unique(tr for bond in everbounded for tr in bond)
x.nb_trajs, len(everbounded_tr), len(everbounded_tr)/float(x.nb_trajs)
sta.plothist(map(len, x.trajs), bins=np.arange(x.size+1))
sta.plothist([len(x.trajs[tr]) for tr in everbounded_tr])
from colloids import statistics as sta
sta.plothist(map(len, x.trajs), bins=np.arange(x.size+1))
sta.plothist([len(x.trajs[tr]) for tr in everbounded_tr])
del _731
sta.__file__
reload(sta)
sta.plothist(map(len, x.trajs), bins=np.arange(x.size+1))
sta.plothist([len(x.trajs[tr]) for tr in everbounded_tr])
everbounded_tr_tot = everbounded_tr[len(x.trajs[tr])==x.size for tr in everbounded_tr]
trajs_tot = [i for i, tr in enumerate(x.trajs) if len(tr)==x.size]
len(trajs_tot), len(everbounded_tr_tot), len(everbounded_tr)/float(len(trajs_tot))
everbounded_tr_tot = everbounded_tr[[len(x.trajs[tr])==x.size for tr in everbounded_tr]]
trajs_tot = [i for i, tr in enumerate(x.trajs) if len(tr)==x.size]
len(trajs_tot), len(everbounded_tr_tot), len(everbounded_tr)/float(len(trajs_tot))
everbounded_tr[0]
everbounded_tr == np.unique(everbounded_tr)
for i,bond in zip(range(10), everbounded):
    print bond
np.unique(tr for i,bond in zip(range(10), everbounded) for tr in bond)
everbounded_tr_tot = everbounded_tr[[len(x.trajs[tr])==x.size for tr in everbounded_tr]]
trajs_tot = [i for i, tr in enumerate(x.trajs) if len(tr)==x.size]
len(trajs_tot), len(everbounded_tr_tot), len(everbounded_tr_tot)/float(len(trajs_tot))
everbounded_tr_tot = [tr for tr everbounded_tr if len(x.trajs[tr])==x.size]
trajs_tot = [i for i, tr in enumerate(x.trajs) if len(tr)==x.size]
len(trajs_tot), len(everbounded_tr_tot), len(everbounded_tr_tot)/float(len(trajs_tot))
everbounded_tr_tot = [tr for tr in everbounded_tr if len(x.trajs[tr])==x.size]
trajs_tot = [i for i, tr in enumerate(x.trajs) if len(tr)==x.size]
len(trajs_tot), len(everbounded_tr_tot), len(everbounded_tr_tot)/float(len(trajs_tot))
get_ipython().magic(u'pinfo %logstart')
get_ipython().magic(u'pwd ()')
get_ipython().magic(u'logstart 2015-autumn.py append')
everbounded_p = np.zeros((x.size, len(everbounded_tr_tot), 3))
for t, name in x.enum():
    pos = np.loadtxt(name, skiprows=2)
    for itr, tr in enumerate(everbounded_tr_tot):
        everbounded_p[t, itr] = pos[x.trajs[tr][t]]
everbounded_p = np.zeros((x.size, len(everbounded_tr_tot), 3))
pro = ProgressBar(x.size)
for t, name in x.enum():
    pos = np.loadtxt(name, skiprows=2)
    for itr, tr in enumerate(everbounded_tr_tot):
        everbounded_p[t, itr] = pos[x.trajs[tr][t]]
    pro.animate(t)
boundnb = np.histogram([tr for bond in everbounded for tr in bond], bins=np.arange(len(x.trajs)))[0]
boundnb.max()
np.argmax(boundnb)
len(everbounded)
tr2everb = -np.ones(len(x.trajs), int)
for itr, tr in enumerate(everbounded_tr_tot):
    tr2everb[tr] = itr
everdists = [
    np.sum(np.diff(everbounded_p[tr2everb[bond]], axis=-1)**2)
    for bond in everbounded
    if tr2everb[bond].min()>-1
    ]    
get_ipython().magic(u'debug')
tr2everb[list(bond)]
tr2everb[bond]
everdists = [
    np.sum(np.diff(everbounded_p[tr2everb[list(bond)]], axis=-1)**2)
    for bond in everbounded
    if tr2everb[list(bond)].min()>-1
    ]    
tr2everb[bond]
tr2everb[list(bond)]
len(tr2everb)
everbounded_p[:,tr2everb[list(bond)]].shape
np.diff(everbounded_p[:,tr2everb[list(bond)]], axis=-2).shape
np.diff(everbounded_p[:,tr2everb[list(bond)]], axis=-2)[:,1].shape
np.diff(everbounded_p[:,tr2everb[list(bond)]], axis=-2)[:,0].shape
everdists = np.array([
    np.sum(np.diff(everbounded_p[:,tr2everb[list(bond)]], axis=-2)[:,0]**2, -1)
    for bond in everbounded
    if tr2everb[list(bond)].min()>-1
    ])
plot(everdists.mean(-1))
plot(everdists.max(-1))
everdists.shape
everdists = np.transpose([
    np.sum(np.diff(everbounded_p[:,tr2everb[list(bond)]], axis=-2)[:,0]**2, -1)
    for bond in everbounded
    if tr2everb[list(bond)].min()>-1
    ])
everdists.shape
everdists = np.array([
    np.sum(np.diff(everbounded_p[:,tr2everb[list(bond)]], axis=-2)[:,0]**2, -1)
    for bond in everbounded
    if tr2everb[list(bond)].min()>-1
    ])
everdists.shape
plot(everdists[0])
everdists = np.sqrt([
    np.sum(np.diff(everbounded_p[:,tr2everb[list(bond)]], axis=-2)[:,0]**2, -1)
    for bond in everbounded
    if tr2everb[list(bond)].min()>-1
    ])
plot(everdists[0])
plot(everdists[1])
plot(everdists[2])
plot(everdists[3])
np.mean(everdists[-1]>1.3*x.radius)
x.radius
1.3*x.radius
2*1.3*x.radius
np.mean(everdists[-1]>2*1.3*x.radius)
plot(np.sqrt(np.mean(everdists[everdists[-1]>2*1.3*x.radius]**2)))
plot(np.sqrt(np.mean(everdists[everdists[-1]>2*1.3*x.radius]**2, 0)))
meandist = np.zeros(x.size)
nb = np.zeros(x.size, int)
for ds in everdists:
    last_att = np.where(ds<2*1.3*x.radius)[-1]
    meandist[:-last_att] += ds[last_att:]**2
    nb[:-last_att] += 1
plot(np.sqrt(meandist / np.maximum(1, nb)))
meandist = np.zeros(x.size)
nb = np.zeros(x.size, int)
for ds in everdists:
    last_att = np.where(ds<2*1.3*x.radius)[0][-1]
    meandist[:-last_att] += ds[last_att:]**2
    nb[:-last_att] += 1
plot(np.sqrt(meandist / np.maximum(1, nb)))
meandist = np.zeros(x.size)
nb = np.zeros(x.size, int)
for ds in everdists:
    last_att = np.where(ds<2*1.3*x.radius)[0][-1]
    meandist[:x.size-last_att] += ds[last_att:]**2
    nb[:x.size-last_att] += 1
plot(np.sqrt(meandist / np.maximum(1, nb)))
plot(np.sqrt(meandist / np.maximum(1, nb)))
xscale('log')
yscale('log')
plot(np.sqrt(meandist / np.maximum(1, nb))/(2*x.rdf_radius()))
xscale('log')
yscale('log')
plot(np.sqrt(meandist / np.maximum(1, nb))/(2*x.rdf_radius()))
tlog = np.logspace(0,3)
plot(tlog, tlog)
plot(tlog, tlog**2)
xscale('log')
yscale('log')
plot(np.sqrt(meandist / np.maximum(1, nb))/(2*x.rdf_radius()))
tlog = np.logspace(0,3)
plot(tlog, tlog)
plot(tlog, tlog**0.5)
xscale('log')
yscale('log')
plot(np.sqrt(meandist / np.maximum(1, nb))/(2*x.rdf_radius()))
tlog = np.logspace(0,3)
plot(tlog, tlog)
plot(tlog, tlog**0.25)
xscale('log')
yscale('log')
plot(np.sqrt(meandist / np.maximum(1, nb))/(2*x.rdf_radius()))
tlog = np.logspace(0,3)
plot(tlog, tlog**0.125)
plot(tlog, tlog**0.25)
xscale('log')
yscale('log')
plot(np.sqrt(meandist / np.maximum(1, nb))/(2*x.rdf_radius()))
tlog = np.logspace(0,3)
plot(tlog, 1.3*tlog**0.125)
plot(tlog, 1.3*tlog**0.25)
xscale('log')
yscale('log')
plot(np.sqrt(meandist / np.maximum(1, nb))/(2*x.rdf_radius()))
tlog = np.logspace(0,3)
plot(tlog, 1.2*tlog**0.125)
plot(tlog, 1.2*tlog**0.25)
xscale('log')
yscale('log')
plot(np.sqrt(meandist / np.maximum(1, nb))/(2*x.rdf_radius()))
tlog = np.logspace(0,3)
plot(tlog, 1.25*tlog**0.125)
plot(tlog, 1.25*tlog**0.1)
xscale('log')
yscale('log')
plot(np.sqrt(meandist / np.maximum(1, nb))/(2*x.rdf_radius()))
tlog = np.logspace(0,3)
#plot(tlog, 1.25*tlog**0.125)
plot(tlog, 1.25*tlog**0.1)
xscale('log')
yscale('log')
plot(meandist / np.maximum(1, nb)/(2*x.rdf_radius()))
tlog = np.logspace(0,3)
#plot(tlog, 1.25*tlog**0.125)
plot(tlog, 1.25**2*tlog**0.2)
xscale('log')
yscale('log')
xlabel('
plot((meandist / np.maximum(1, nb))/(2*x.rdf_radius())**2)
tlog = np.logspace(0,3)
#plot(tlog, 1.25*tlog**0.125)
plot(tlog, 1.25**2*tlog**0.2)
xscale('log')
yscale('log')
xlabel('
plot((meandist / np.maximum(1, nb))/(2*x.rdf_radius())**2)
tlog = np.logspace(0,3)
#plot(tlog, 1.25*tlog**0.125)
plot(tlog, 1.25**2*tlog**0.2)
xscale('log')
yscale('log')
xlabel('t')
ylabel(r'$\Delta r^2/\sigma^2$')
plot((meandist / np.maximum(1, nb))/(2*x.rdf_radius())**2)
tlog = np.logspace(0,3)
#plot(tlog, 1.25*tlog**0.125)
plot(tlog, 1.2**2*tlog**0.2)
xscale('log')
yscale('log')
xlabel('t')
ylabel(r'$\Delta r^2/\sigma^2$')
plot((meandist / np.maximum(1, nb))/(2*x.rdf_radius())**2)
tlog = np.logspace(0,3)
#plot(tlog, 1.25*tlog**0.125)
plot(tlog, 1.22**2*tlog**0.2)
xscale('log')
yscale('log')
xlabel('t')
ylabel(r'$\Delta r^2/\sigma^2$')
def MSDpostbreak(x, L=3):
    """look for all the trajectories that are ever bounded together and unbounded at the next time step
    with a post breaking distence larger than L"""
    pro = ProgressBar(x.size)
    everbounded = set()
    for t, name in x.enum(ext='bonds'):
        try:
            bonds1 = np.atleast_2d(np.loadtxt(name, dtype=int))
        except UserWarning:
            bonds1 = np.zeros([0,2], int)
        p2tr1 = np.loadtxt(x.get_format_string(ext='p2tr')%t, dtype=int)
        if t>0:
            everbounded |= set(
                tuple(sorted((path[0], path[-1]))) 
                for path in broken_bonds_path(bonds0, bonds1, p2tr0, p2tr1)
                if len(path)>L+1
                )
        bonds0 = bonds1
        p2tr0 = p2tr1
        pro.animate(t)
    #single trajectories that are bounded
    everbounded_tr = np.unique(tr for bond in everbounded for tr in bond)
    #filter out the trajectories that do not span the total length
    everbounded_tr_tot = [tr for tr in everbounded_tr if len(x.trajs[tr])==x.size]
    #index them
    tr2everb = -np.ones(len(x.trajs), int)
    for itr, tr in enumerate(everbounded_tr_tot):
        tr2everb[tr] = itr
    #load positions
    everbounded_p = np.zeros((x.size, len(everbounded_tr_tot), 3))
    pro = ProgressBar(x.size)
    for t, name in x.enum():
        pos = np.loadtxt(name, skiprows=2)
        for itr, tr in enumerate(everbounded_tr_tot):
            everbounded_p[t, itr] = pos[x.trajs[tr][t]]
        pro.animate(t)
    #distances between the couples of trajectories
    everdists = np.sqrt([
        np.sum(np.diff(everbounded_p[:,tr2everb[list(bond)]], axis=-2)[:,0]**2, -1)
        for bond in everbounded
        if tr2everb[list(bond)].min()>-1
        ])
    #mean square distance from the last time bounded
    meandist = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for ds in everdists:
        last_att = np.where(ds<2*1.3*x.radius)[0][-1]
        meandist[:x.size-last_att] += ds[last_att:]**2
        nb[:x.size-last_att] += 1
    return (meandist / np.maximum(1, nb))/(2*x.rdf_radius())**2
msds = [MSDpostbreak(x, L) for L in range(2,10)]
for L, msd in zip(range(2,10), msds):
    plot(msd, label=L)
for L, msd in zip(range(2,10), msds):
    plot(msd, label=L, cm.autumn(0.1*L9))
xscale('log')
yscale('log')
xlabel('t')
ylabel(r'$\Delta r^2/\sigma^2$')
for L, msd in zip(range(2,10), msds):
    plot(msd, label=L, color=cm.autumn(0.1*L9))
xscale('log')
yscale('log')
xlabel('t')
ylabel(r'$\Delta r^2/\sigma^2$')
for L, msd in zip(range(2,10), msds):
    plot(msd, label=L, color=cm.autumn(0.1*L))
xscale('log')
yscale('log')
xlabel('t')
ylabel(r'$\Delta r^2/\sigma^2$')
for L, msd in zip(range(2,10), msds):
    plot(msd, label=L, color=cm.autumn(0.1*L))
tlog = np.logspace(0,3)
plot(tlog, 1.22**2*tlog**0.2, color='k')
xscale('log')
yscale('log')
xlabel('t')
ylabel(r'$\Delta r^2/\sigma^2$')
for L, msd in zip(range(2,10), msds):
    plot(msd, label=L, color=cm.autumn(0.1*L))
tlog = np.logspace(0,3)
plot(tlog, 1.22**2*tlog**0.2, color='k')
plot(tlog, 1.3**2*tlog**0.2, color='k')
xscale('log')
yscale('log')
xlabel('t')
ylabel(r'$\Delta r^2/\sigma^2$')
for L, msd in zip(range(2,10), msds):
    plot(msd, label=L, color=cm.autumn(0.1*L))
tlog = np.logspace(0,3)
plot(tlog, 1.22**2*tlog**0.2, color='k')
plot(tlog, 1.3**2*tlog**0.5, color='k')
xscale('log')
yscale('log')
xlabel('t')
ylabel(r'$\Delta r^2/\sigma^2$')
for L, msd in zip(range(2,10), msds):
    plot(msd, label=L, color=cm.autumn(0.1*L))
tlog = np.logspace(0,3)
plot(tlog, 1.22**2*tlog**0.2, color='k')
plot(tlog, 1.3**2*tlog**0.4, color='k')
xscale('log')
yscale('log')
xlabel('t')
ylabel(r'$\Delta r^2/\sigma^2$')
for L, msd in zip(range(2,10), msds):
    plot(msd, label=L, color=cm.autumn(0.1*L))
tlog = np.logspace(0,3)
plot(tlog, 1.22**2*tlog**0.2, color='k')
plot(tlog, 1.3**2*tlog**0.3, color='k')
xscale('log')
yscale('log')
xlabel('t')
ylabel(r'$\Delta r^2/\sigma^2$')
for L, msd in zip(range(2,10), msds):
    plot(msd, label=L, color=cm.autumn(0.1*L))
tlog = np.logspace(0,3)
plot(tlog, 1.22**2*tlog**0.2, color='k')
plot(tlog, 1.22**2*tlog**0.35, color='k')
xscale('log')
yscale('log')
xlabel('t')
ylabel(r'$\Delta r^2/\sigma^2$')
for L, msd in zip(range(2,10), msds):
    plot(msd, label=L, color=cm.autumn(0.1*L))
tlog = np.logspace(0,3)
plot(tlog, 1.22**2*tlog**0.2, color='k')
plot(tlog, 1.22**2*tlog**0.36, color='k')
xscale('log')
yscale('log')
xlabel('t')
ylabel(r'$\Delta r^2/\sigma^2$')
plot((meandist / np.maximum(1, nb))/(2*x.rdf_radius())**2)
for L, msd in zip(range(2,10), msds):
    plot(msd, label=L, color=cm.autumn(0.1*L))
tlog = np.logspace(0,3)
plot(tlog, 1.22**2*tlog**0.2, color='k')
plot(tlog, 1.22**2*tlog**0.36, color='k')
xscale('log')
yscale('log')
xlabel('t')
ylabel(r'$\Delta r^2/\sigma^2$')
plot((meandist / np.maximum(1, nb))/(2*x.rdf_radius())**2)
for L, msd in zip(range(2,10), msds):
    plot(msd, label=L, color=cm.autumn(0.1*L))
tlog = np.logspace(0,3)
plot(tlog, 1.22**2*tlog**0.2, color='k')
plot(tlog, 1.22**2*tlog**0.38, color='k')
xscale('log')
yscale('log')
xlabel('t')
ylabel(r'$\Delta r^2/\sigma^2$')
np.save(os.path.join(x.path, 'broken_msd.npy'), np.column_stack(((meandist / np.maximum(1, nb))/(2*x.rdf_radius())**2, msds)))
meandist.shape, msds.shape
np.save(os.path.join(x.path, 'broken_msd.npy'), np.column_stack([(meandist / np.maximum(1, nb))/(2*x.rdf_radius())**2] + msds))
plot((meandist / np.maximum(1, nb))/(2*x.rdf_radius())**2)
for L, msd in zip(range(2,10), msds):
    plot(msd, label=L, color=cm.autumn(0.1*L))
tlog = np.logspace(0,3)
plot(tlog, 1.22**2*tlog**0.2, color='k')
plot(tlog, 1.22**2*tlog**0.38, color='k')
xscale('log')
yscale('log')
ylim(1,10)
xlabel('t')
ylabel(r'$\Delta r^2/\sigma^2$')
plot((meandist / np.maximum(1, nb))/(2*x.rdf_radius())**2)
for L, msd in zip(range(2,10), msds):
    plot(msd, label=L, color=cm.autumn(0.1*L))
tlog = np.logspace(0,3)
plot(tlog, 1.22**2*tlog**0.2, color='k')
plot(tlog, 1.22**2*tlog**0.38, color='k')
xscale('log')
yscale('log')
ylim(1,10)
xlim(1,100)
xlabel('t')
ylabel(r'$\Delta r^2/\sigma^2$')
plot((meandist / np.maximum(1, nb))/(2*x.rdf_radius())**2)
for L, msd in zip(range(2,10), msds):
    plot(msd, label=L, color=cm.autumn(0.1*L))
tlog = np.logspace(0,3)
plot(tlog, 1.22**2*tlog**0.2, color='k')
plot(tlog, 1.22**2*tlog**0.38, color='k')
xscale('log')
yscale('log')
ylim(1,10)
xlim(1,x.size-68)
xlabel('t')
ylabel(r'$\Delta r^2/\sigma^2$')
plot((meandist / np.maximum(1, nb))/(2*x.rdf_radius())**2)
for L, msd in zip(range(2,10), msds):
    plot(msd, label=L, color=cm.autumn(0.1*L))
tlog = np.logspace(0,3)
plot(tlog, 1.22**2*tlog**0.2, color='k')
plot(tlog, 1.22**2*tlog**0.38, color='k')
xscale('log')
yscale('log')
ylim(1,10)
xlim(1,x.size-69)
xlabel('t')
ylabel(r'$\Delta r^2/\sigma^2$')
plot((meandist / np.maximum(1, nb))/(2*x.rdf_radius())**2)
for L, msd in zip(range(2,10), msds):
    plot(msd, label=L, color=cm.autumn(0.1*L))
tlog = np.logspace(0,3)
plot(tlog, 1.22**2*tlog**0.2, color='k')
plot(tlog, 1.22**2*tlog**0.38, color='k')
xscale('log')
yscale('log')
ylim(1,10)
xlim(1,x.size-70)
xlabel('t')
ylabel(r'$\Delta r^2/\sigma^2$')
plot((meandist / np.maximum(1, nb))/(2*x.rdf_radius())**2)
for L, msd in zip(range(2,10), msds):
    plot(msd, label=L, color=cm.autumn(0.1*L))
tlog = np.logspace(0,3)
plot(tlog, 1.22**2*tlog**0.2, color='k')
plot(tlog, 1.22**2*tlog**0.38, color='k')
xscale('log')
yscale('log')
ylim(1,10)
xlim(1,x.size-70)
xlabel(r't/\tau_B')
ylabel(r'$\Delta r^2/\sigma^2$')
savefig(os.path.join(x.path, 'broken_msd.pdf'))
plot((meandist / np.maximum(1, nb))/(2*x.rdf_radius())**2)
for L, msd in zip(range(2,10), msds):
    plot(msd, label=L, color=cm.autumn(0.1*L))
tlog = np.logspace(0,3)
plot(tlog, 1.22**2*tlog**0.2, color='k')
plot(tlog, 1.22**2*tlog**0.38, color='k')
xscale('log')
yscale('log')
ylim(1,10)
xlim(1,x.size-70)
xlabel(r'$t/\tau_B$')
ylabel(r'$\Delta r^2/\sigma^2$')
savefig(os.path.join(x.path, 'broken_msd.pdf'))
plot((meandist / np.maximum(1, nb))/(2*x.rdf_radius())**2)
for L, msd in zip(range(2,10), msds):
    plot(msd, label=L, color=cm.autumn(0.1*L))
tlog = np.logspace(0,3)
plot(tlog, 1.22**2*tlog**0.2, color='k')
plot(tlog, 1.22**2*tlog**0.38, color='k')
xscale('log')
yscale('log')
ylim(1,10)
xlim(1,x.size-70)
xlabel(r'$\Delta t/\tau_B$')
ylabel(r'$\Delta r^2/\sigma^2$')
savefig(os.path.join(x.path, 'broken_msd.pdf'))
os.path.join(x.path, 'broken_msd.pdf')
nbb = np.zeros(x.size, int)
for ds in everdists:
    last_att = np.where(ds<2*1.3*x.radius)[0][-1]
    nbb[last_att] += 1
plot(nbb[t0:])
xscale('log')
yscale('log')
msds = [MSDpostbreak(x, L) for L in [0]+range(2,10)]
#plot((meandist / np.maximum(1, nb))/(2*x.rdf_radius())**2)
for L, msd in zip(range(0,2,10), msds):
    plot(msd, label=L, color=cm.autumn(0.1*L))
tlog = np.logspace(0,3)
plot(tlog, 1.22**2*tlog**0.2, color='k')
plot(tlog, 1.22**2*tlog**0.38, color='k')
xscale('log')
yscale('log')
ylim(1,10)
xlim(1,x.size-70)
xlabel(r'$\Delta t/\tau_B$')
ylabel(r'$\Delta r^2/\sigma^2$')
savefig(os.path.join(x.path, 'broken_msd.pdf'))
#plot((meandist / np.maximum(1, nb))/(2*x.rdf_radius())**2)
for L, msd in zip([0]+range(2,10), msds):
    plot(msd, label=L, color=cm.autumn(0.1*L))
tlog = np.logspace(0,3)
plot(tlog, 1.22**2*tlog**0.2, color='k')
plot(tlog, 1.22**2*tlog**0.38, color='k')
xscale('log')
yscale('log')
ylim(1,10)
xlim(1,x.size-70)
xlabel(r'$\Delta t/\tau_B$')
ylabel(r'$\Delta r^2/\sigma^2$')
savefig(os.path.join(x.path, 'broken_msd.pdf'))
def neverreform(x, L=3, mintrajlength=None):
    """All the broken bonds (in term of trajectories) with a post breaking distence larger than L, assiciated with and the last time they break"""
    pro = ProgressBar(x.size)
    result = dict()
    for t, name in x.enum(ext='bonds'):
        try:
            bonds1 = np.atleast_2d(np.loadtxt(name, dtype=int))
        except UserWarning:
            bonds1 = np.zeros([0,2], int)
        p2tr1 = np.loadtxt(x.get_format_string(ext='p2tr')%t, dtype=int)
        if t>0:
            for path in broken_bonds_path(bonds0, bonds1, p2tr0, p2tr1):
                if len(path)<L:
                    continue
                a, b = sorted([path[0], path[-1]])
                if mintrajlength is not None and min(len(x.trajs[u]) for u in [a,b]) <= mintrajlength:
                    continue
                result[(a,b)] = t #create or replace if was already broken before
        bonds0 = bonds1
        p2tr0 = p2tr1
        pro.animate(t)
    return result
nev = neverreform(x)
plot(nbb[t0:])
plot(np.histogram(nev.itervalues(), np.arange(x.size+1))[0])
xscale('log')
yscale('log')
plot(nbb[t0:])
plot(np.histogram(nev.values(), np.arange(t0, x.size+1))[0])
xscale('log')
yscale('log')
plot(nbb[t0:])
plot(np.histogram(nev.values(), np.arange(t0, x.size+1))[0])
plot(np.histogram([
    t for (a,b) in nev.iteritems()
    if min(len(x.trajs[u]) for u in [a,b]) == x.size
    ], np.arange(t0, x.size+1))[0])
xscale('log')
yscale('log')
plot(nbb[t0:])
plot(np.histogram(nev.values(), np.arange(t0, x.size+1))[0])
plot(np.histogram([
    t for (a,b),t in nev.iteritems()
    if min(len(x.trajs[u]) for u in [a,b]) == x.size
    ], np.arange(t0, x.size+1))[0])
xscale('log')
yscale('log')
plot(nbb[t0:])
plot(np.histogram(nev.values(), np.arange(t0+1, x.size+1))[0])
plot(np.histogram([
    t for (a,b),t in nev.iteritems()
    if min(len(x.trajs[u]) for u in [a,b]) == x.size
    ], np.arange(t0+1, x.size+1))[0])
xscale('log')
yscale('log')
plot(nbb[t0:])
plot(np.histogram(nev.values(), np.arange(t0+1, x.size+1))[0])
plot(np.histogram([
    t for (a,b),t in nev.iteritems()
    if min(len(x.trajs[u]) for u in [a,b]) == x.size
    ], np.arange(t0+1, x.size+1))[0])
plot([10,100], np.power([10,100], 1.5))
xscale('log')
yscale('log')
plot(nbb[t0:])
plot(np.histogram(nev.values(), np.arange(t0+1, x.size+1))[0])
plot(np.histogram([
    t for (a,b),t in nev.iteritems()
    if min(len(x.trajs[u]) for u in [a,b]) == x.size
    ], np.arange(t0+1, x.size+1))[0])
plot([10,100], np.power([10,100], -1.5))
xscale('log')
yscale('log')
plot(nbb[t0:])
plot(np.histogram(nev.values(), np.arange(t0+1, x.size+1))[0])
plot(np.histogram([
    t for (a,b),t in nev.iteritems()
    if min(len(x.trajs[u]) for u in [a,b]) == x.size
    ], np.arange(t0+1, x.size+1))[0])
plot([10,100], 100*np.power([10,100], -1.5))
xscale('log')
yscale('log')
plot(nbb[t0:])
plot(np.histogram(nev.values(), np.arange(t0+1, x.size+1))[0])
plot(np.histogram([
    t for (a,b),t in nev.iteritems()
    if min(len(x.trajs[u]) for u in [a,b]) == x.size
    ], np.arange(t0+1, x.size+1))[0])
plot([10,100], 10e4*np.power([10,100], -1.5))
xscale('log')
yscale('log')
plot(nbb[t0:])
plot(np.histogram(nev.values(), np.arange(t0+1, x.size+1))[0])
plot(np.histogram([
    t for (a,b),t in nev.iteritems()
    if min(len(x.trajs[u]) for u in [a,b]) == x.size
    ], np.arange(t0+1, x.size+1))[0])
plot([10,100], 5e4*np.power([10,100], -1.5))
xscale('log')
yscale('log')
plot(nbb[t0:])
plot(np.histogram(nev.values(), np.arange(t0+1, x.size+1))[0])
plot(np.histogram([
    t for (a,b),t in nev.iteritems()
    if min(len(x.trajs[u]) for u in [a,b]) == x.size
    ], np.arange(t0+1, x.size+1))[0])
plot([10,100], 2e4*np.power([10,100], -1.5))
xscale('log')
yscale('log')
plot(nbb[t0:])
plot(np.histogram(nev.values(), np.arange(t0+1, x.size+1))[0])
plot(np.histogram([
    t for (a,b),t in nev.iteritems()
    if min(len(x.trajs[u]) for u in [a,b]) == x.size
    ], np.arange(t0+1, x.size+1))[0])
plot([10,100], 2e4*np.power([10,100], -1.5))
xscale('log')
yscale('log')
plot(nbb[t0:]*1./x.get_nb()[t0:-1])
plot(np.histogram(nev.values(), np.arange(t0+1, x.size+1))[0]*1./x.get_nb()[t0:-1])
plot(np.histogram(
    [
        t for (a,b),t in nev.iteritems()
        if min(len(x.trajs[u]) for u in [a,b]) == x.size
        ], 
    np.arange(t0+1, x.size+1))[0]*1./x.get_nb()[t0:-1])

maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
for L in [5,13, 23]:
    rates = []
    for head,t0 in zip(['1_percolation/1340_percolation', '2_ageing/1415_ageing'], [68, 0]):
        x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/%s.traj'%(head))
        lengths = [
            np.loadtxt(name, int) 
            for t, name in x.enum('_broken', 'length') 
            if t>0]
        maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
        his = np.array([
            np.histogram(l, np.arange(-1, maxlength+1))[0]
            for l in lengths])
        rates.append(np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )))
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], label=lab)[0]



plot([10,100], 2e4*np.power([10,100], -1.5))
xscale('log')
yscale('log')
plot(nbb[t0:]*1./x.get_nb()[t0:])
plot(np.histogram(nev.values(), np.arange(t0+1, x.size+1))[0]*1./x.get_nb()[t0:-1])
plot(np.histogram(
    [
        t for (a,b),t in nev.iteritems()
        if min(len(x.trajs[u]) for u in [a,b]) == x.size
        ], 
    np.arange(t0+1, x.size+1))[0]*1./x.get_nb()[t0:-1])

maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
for L in [5,13, 23]:
    rates = []
    for head,t0 in zip(['1_percolation/1340_percolation', '2_ageing/1415_ageing'], [68, 0]):
        x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/%s.traj'%(head))
        lengths = [
            np.loadtxt(name, int) 
            for t, name in x.enum('_broken', 'length') 
            if t>0]
        maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
        his = np.array([
            np.histogram(l, np.arange(-1, maxlength+1))[0]
            for l in lengths])
        rates.append(np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )))
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], label=lab)[0]



plot([10,100], 2e4*np.power([10,100], -1.5))
xscale('log')
yscale('log')
plot(nbb[t0:]*1./x.get_nb()[t0:])
plot(np.histogram(nev.values(), np.arange(t0+1, x.size+1))[0]*1./x.get_nb()[t0:-1])
plot(np.histogram(
    [
        t for (a,b),t in nev.iteritems()
        if min(len(x.trajs[u]) for u in [a,b]) == x.size
        ], 
    np.arange(t0+1, x.size+1))[0]*1./x.get_nb()[t0:-1])

for L in [5,13, 23]:
    rates = []
    for head,t0 in zip(['1_percolation/1340_percolation', '2_ageing/1415_ageing'], [68, 0]):
        x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/%s.traj'%(head))
        lengths = [
            np.loadtxt(name, int) 
            for t, name in x.enum('_broken', 'length') 
            if t>0]
        maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
        his = np.array([
            np.histogram(l, np.arange(-1, maxlength+1))[0]
            for l in lengths])
        rates.append(np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )))
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], label=lab)[0]



plot([10,100], 2e4*np.power([10,100], -1.5))
xscale('log')
yscale('log')
plot(nbb[t0:]*1./x.get_nb()[t0:])
plot(np.histogram(nev.values(), np.arange(t0+1, x.size+1))[0]*1./x.get_nb()[t0:-1])
plot(np.histogram(
    [
        t for (a,b),t in nev.iteritems()
        if min(len(x.trajs[u]) for u in [a,b]) == x.size
        ], 
    np.arange(t0+1, x.size+1))[0]*1./x.get_nb()[t0:-1])

for L in [5,13, 23]:
    rates = []
    for head,t0 in zip(['1_percolation/1340_percolation', '2_ageing/1415_ageing'], [68, 0]):
        x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/%s.traj'%(head))
        lengths = [
            np.loadtxt(name, int) 
            for t, name in x.enum('_broken', 'length') 
            if t>0]
        maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
        his = np.array([
            np.histogram(l, np.arange(-1, maxlength+1))[0]
            for l in lengths])
        rates.append(np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )))
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], label=lab)[0]



plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
xscale('log')
yscale('log')
x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
plot(nbb[t0:]*1./x.get_nb()[t0:])
plot(np.histogram(nev.values(), np.arange(t0+1, x.size+1))[0]*1./x.get_nb()[t0:-1])
plot(np.histogram(
    [
        t for (a,b),t in nev.iteritems()
        if min(len(x.trajs[u]) for u in [a,b]) == x.size
        ], 
    np.arange(t0+1, x.size+1))[0]*1./x.get_nb()[t0:-1])

for L in [5,13, 23]:
    rates = []
    for head,t0 in zip(['1_percolation/1340_percolation', '2_ageing/1415_ageing'], [68, 0]):
        x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/%s.traj'%(head))
        lengths = [
            np.loadtxt(name, int) 
            for t, name in x.enum('_broken', 'length') 
            if t>0]
        maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
        his = np.array([
            np.histogram(l, np.arange(-1, maxlength+1))[0]
            for l in lengths])
        rates.append(np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )))
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], label=lab)[0]



plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
xscale('log')
yscale('log')
x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
t0 = 68
plot(nbb[t0:]*1./x.get_nb()[t0:])
plot(np.histogram(nev.values(), np.arange(t0+1, x.size+1))[0]*1./x.get_nb()[t0:-1])
plot(np.histogram(
    [
        t for (a,b),t in nev.iteritems()
        if min(len(x.trajs[u]) for u in [a,b]) == x.size
        ], 
    np.arange(t0+1, x.size+1))[0]*1./x.get_nb()[t0:-1])

for L in [5,13, 23]:
    rates = []
    for head,t0 in zip(['1_percolation/1340_percolation', '2_ageing/1415_ageing'], [68, 0]):
        x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/%s.traj'%(head))
        lengths = [
            np.loadtxt(name, int) 
            for t, name in x.enum('_broken', 'length') 
            if t>0]
        maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
        his = np.array([
            np.histogram(l, np.arange(-1, maxlength+1))[0]
            for l in lengths])
        rates.append(np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )))
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], label=lab)[0]



plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
xscale('log')
yscale('log')
x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
t0 = 68
plot(nbb[t0+1:]*1./x.get_nb()[t0:-1])
plot(np.histogram(nev.values(), np.arange(t0+1, x.size+1))[0]*1./x.get_nb()[t0:-1])
plot(np.histogram(
    [
        t for (a,b),t in nev.iteritems()
        if min(len(x.trajs[u]) for u in [a,b]) == x.size
        ], 
    np.arange(t0+1, x.size+1))[0]*1./x.get_nb()[t0:-1])

for L in [5,13, 23]:
    rates = []
    for head,t0 in zip(['1_percolation/1340_percolation', '2_ageing/1415_ageing'], [68, 0]):
        x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/%s.traj'%(head))
        lengths = [
            np.loadtxt(name, int) 
            for t, name in x.enum('_broken', 'length') 
            if t>0]
        maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
        his = np.array([
            np.histogram(l, np.arange(-1, maxlength+1))[0]
            for l in lengths])
        rates.append(np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )))
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], label=lab)[0]



plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
xscale('log')
yscale('log')
x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
nev = neverreform(x,4)
x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/2_ageing/1415_ageing.traj')
nev2 = neverreform(x,4)
x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
t0 = 68
plot(nbb[t0+1:]*1./x.get_nb()[t0:-1])
plot(np.histogram(nev.values(), np.arange(t0+1, x.size+1))[0]*1./x.get_nb()[t0:-1])
plot(np.histogram(
    [
        t for (a,b),t in nev.iteritems()
        if min(len(x.trajs[u]) for u in [a,b]) == x.size
        ], 
    np.arange(t0+1, x.size+1))[0]*1./x.get_nb()[t0:-1])

for L in [5,13, 23]:
    rates = []
    for head,t0 in zip(['1_percolation/1340_percolation', '2_ageing/1415_ageing'], [68, 0]):
        x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/%s.traj'%(head))
        lengths = [
            np.loadtxt(name, int) 
            for t, name in x.enum('_broken', 'length') 
            if t>0]
        maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
        his = np.array([
            np.histogram(l, np.arange(-1, maxlength+1))[0]
            for l in lengths])
        rates.append(np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )))
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], label=lab)[0]



plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
xscale('log')
yscale('log')
r_nev = np.vstack([
    np.column_stack((
        np.arange(t0+1, x.size)*Dt+t00,
        np.histogram(n.values(), np.arange(t0+1, x.size+1))[0]*1./x.get_nb()[t0:-1]
        ))
    for x, n, t0, Dt, t00 in zip([xs[0], xs[-3]], [nev, nev2], [68,0], [1,3], [0, 205-68])
    ])
x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
t0 = 68
plot(nbb[t0+1:]*1./x.get_nb()[t0:-1])
#plot(np.histogram(nev.values(), np.arange(t0+1, x.size+1))[0]*1./x.get_nb()[t0:-1])
plot(r_nev[:,0], r_nev[:,1])
plot(np.histogram(
    [
        t for (a,b),t in nev.iteritems()
        if min(len(x.trajs[u]) for u in [a,b]) == x.size
        ], 
    np.arange(t0+1, x.size+1))[0]*1./x.get_nb()[t0:-1])

for L in [5,13, 23]:
    rates = []
    for head,t0 in zip(['1_percolation/1340_percolation', '2_ageing/1415_ageing'], [68, 0]):
        x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/%s.traj'%(head))
        lengths = [
            np.loadtxt(name, int) 
            for t, name in x.enum('_broken', 'length') 
            if t>0]
        maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
        his = np.array([
            np.histogram(l, np.arange(-1, maxlength+1))[0]
            for l in lengths])
        rates.append(np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )))
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], label=lab)[0]



plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
xscale('log')
yscale('log')
r_nev = np.vstack([
    np.column_stack((
        np.arange(x.size)*Dt+t00,
        np.histogram(n.values(), np.arange(x.size+1))[0]*1./x.get_nb()[:-1]
        ))
    for x, n, Dt, t00 in zip([xs[0], xs[-3]], [nev, nev2], [1,3], [0, 205-68])
    ])
r_nev = np.vstack([
    np.column_stack((
        np.arange(x.size)*Dt+t00,
        np.histogram(n.values(), np.arange(x.size))[0]*1./x.get_nb()[:-1]
        ))
    for x, n, Dt, t00 in zip([xs[0], xs[-3]], [nev, nev2], [1,3], [0, 205-68])
    ])
r_nev = np.vstack([
    np.column_stack((
        np.arange(x.size-1)*Dt+t00,
        np.histogram(n.values(), np.arange(x.size))[0]*1./x.get_nb()[:-1]
        ))
    for x, n, Dt, t00 in zip([xs[0], xs[-3]], [nev, nev2], [1,3], [0, 205-68])
    ])
x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
t0 = 68
plot(nbb[t0+1:]*1./x.get_nb()[t0:-1])
#plot(np.histogram(nev.values(), np.arange(t0+1, x.size+1))[0]*1./x.get_nb()[t0:-1])
plot(r_nev[:,0], r_nev[:,1])
plot(np.histogram(
    [
        t for (a,b),t in nev.iteritems()
        if min(len(x.trajs[u]) for u in [a,b]) == x.size
        ], 
    np.arange(t0+1, x.size+1))[0]*1./x.get_nb()[t0:-1])

for L in [5,13, 23]:
    rates = []
    for head,t0 in zip(['1_percolation/1340_percolation', '2_ageing/1415_ageing'], [68, 0]):
        x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/%s.traj'%(head))
        lengths = [
            np.loadtxt(name, int) 
            for t, name in x.enum('_broken', 'length') 
            if t>0]
        maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
        his = np.array([
            np.histogram(l, np.arange(-1, maxlength+1))[0]
            for l in lengths])
        rates.append(np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )))
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], label=lab)[0]



plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
xscale('log')
yscale('log')
x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
t0 = 68
plot(nbb[t0+1:]*1./x.get_nb()[t0:-1])
#plot(np.histogram(nev.values(), np.arange(t0+1, x.size+1))[0]*1./x.get_nb()[t0:-1])
plot(r_nev[t0:,0], r_nev[:,1])
plot(np.histogram(
    [
        t for (a,b),t in nev.iteritems()
        if min(len(x.trajs[u]) for u in [a,b]) == x.size
        ], 
    np.arange(t0+1, x.size+1))[0]*1./x.get_nb()[t0:-1])

for L in [5,13, 23]:
    rates = []
    for head,t0 in zip(['1_percolation/1340_percolation', '2_ageing/1415_ageing'], [68, 0]):
        x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/%s.traj'%(head))
        lengths = [
            np.loadtxt(name, int) 
            for t, name in x.enum('_broken', 'length') 
            if t>0]
        maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
        his = np.array([
            np.histogram(l, np.arange(-1, maxlength+1))[0]
            for l in lengths])
        rates.append(np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )))
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], label=lab)[0]



plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
xscale('log')
yscale('log')
x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
t0 = 68
plot(nbb[t0+1:]*1./x.get_nb()[t0:-1])
#plot(np.histogram(nev.values(), np.arange(t0+1, x.size+1))[0]*1./x.get_nb()[t0:-1])
plot(r_nev[t0:,0], r_nev[t0:,1])
plot(np.histogram(
    [
        t for (a,b),t in nev.iteritems()
        if min(len(x.trajs[u]) for u in [a,b]) == x.size
        ], 
    np.arange(t0+1, x.size+1))[0]*1./x.get_nb()[t0:-1])

for L in [5,13, 23]:
    rates = []
    for head,t0 in zip(['1_percolation/1340_percolation', '2_ageing/1415_ageing'], [68, 0]):
        x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/%s.traj'%(head))
        lengths = [
            np.loadtxt(name, int) 
            for t, name in x.enum('_broken', 'length') 
            if t>0]
        maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
        his = np.array([
            np.histogram(l, np.arange(-1, maxlength+1))[0]
            for l in lengths])
        rates.append(np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )))
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], label=lab)[0]



plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
xscale('log')
yscale('log')
r_nev = np.vstack([
    np.column_stack((
        np.arange(x.size-1)*Dt+t00,
        np.histogram(n.values(), np.arange(x.size))[0]*1./x.get_nb()[:-1]
        ))
    for x, n, Dt, t00 in zip([xs[0], xs[-3]], [nev, nev2], [1,3], [0, 205])
    ])
x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
t0 = 68
plot(nbb[t0+1:]*1./x.get_nb()[t0:-1])
#plot(np.histogram(nev.values(), np.arange(t0+1, x.size+1))[0]*1./x.get_nb()[t0:-1])
plot(r_nev[t0:,0]-t0, r_nev[t0:,1])
plot(np.histogram(
    [
        t for (a,b),t in nev.iteritems()
        if min(len(x.trajs[u]) for u in [a,b]) == x.size
        ], 
    np.arange(t0+1, x.size+1))[0]*1./x.get_nb()[t0:-1])

for L in [5,13, 23]:
    rates = []
    for head,t0 in zip(['1_percolation/1340_percolation', '2_ageing/1415_ageing'], [68, 0]):
        x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/%s.traj'%(head))
        lengths = [
            np.loadtxt(name, int) 
            for t, name in x.enum('_broken', 'length') 
            if t>0]
        maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
        his = np.array([
            np.histogram(l, np.arange(-1, maxlength+1))[0]
            for l in lengths])
        rates.append(np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )))
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], label=lab)[0]



plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
xscale('log')
yscale('log')
x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
t0 = 68
plot(nbb[t0+1:]*1./x.get_nb()[t0:-1])
#plot(np.histogram(nev.values(), np.arange(t0+1, x.size+1))[0]*1./x.get_nb()[t0:-1])
plot(r_nev[t0:,0]-t0-1, r_nev[t0:,1])
plot(np.histogram(
    [
        t for (a,b),t in nev.iteritems()
        if min(len(x.trajs[u]) for u in [a,b]) == x.size
        ], 
    np.arange(t0+1, x.size+1))[0]*1./x.get_nb()[t0:-1])

for L in [5,13, 23]:
    rates = []
    for head,t0 in zip(['1_percolation/1340_percolation', '2_ageing/1415_ageing'], [68, 0]):
        x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/%s.traj'%(head))
        lengths = [
            np.loadtxt(name, int) 
            for t, name in x.enum('_broken', 'length') 
            if t>0]
        maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
        his = np.array([
            np.histogram(l, np.arange(-1, maxlength+1))[0]
            for l in lengths])
        rates.append(np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )))
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], label=lab)[0]



plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
xscale('log')
yscale('log')
r_revspan = np.vstack([
    np.column_stack((
        np.arange(x.size-1)*Dt+t00,
        np.histogram(
            [
                t for (a,b),t in n.iteritems()
                if min(len(x.trajs[u]) for u in [a,b]) == x.size
                ], 
            np.arange(x.size)
            )[0]*1./x.get_nb()[:-1]
        ))
    for x, n, Dt, t00 in zip([xs[0], xs[-3]], [nev, nev2], [1,3], [0, 205])
    ])
x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
t0 = 68
plot(nbb[t0+1:]*1./x.get_nb()[t0:-1])
#plot(np.histogram(nev.values(), np.arange(t0+1, x.size+1))[0]*1./x.get_nb()[t0:-1])
plot(r_nev[t0:,0]-t0-1, r_nev[t0:,1])
#plot(np.histogram(
 #   [
  #      t for (a,b),t in nev.iteritems()
   #     if min(len(x.trajs[u]) for u in [a,b]) == x.size
    #    ], 
    #np.arange(t0+1, x.size+1))[0]*1./x.get_nb()[t0:-1])
plot(r_nevspan[t0:,0]-t0-1, r_nevspan[t0:,1])

for L in [5,13, 23]:
    rates = []
    for head,t0 in zip(['1_percolation/1340_percolation', '2_ageing/1415_ageing'], [68, 0]):
        x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/%s.traj'%(head))
        lengths = [
            np.loadtxt(name, int) 
            for t, name in x.enum('_broken', 'length') 
            if t>0]
        maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
        his = np.array([
            np.histogram(l, np.arange(-1, maxlength+1))[0]
            for l in lengths])
        rates.append(np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )))
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], label=lab)[0]



plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
xscale('log')
yscale('log')
r_nevspan = np.vstack([
    np.column_stack((
        np.arange(x.size-1)*Dt+t00,
        np.histogram(
            [
                t for (a,b),t in n.iteritems()
                if min(len(x.trajs[u]) for u in [a,b]) == x.size
                ], 
            np.arange(x.size)
            )[0]*1./x.get_nb()[:-1]
        ))
    for x, n, Dt, t00 in zip([xs[0], xs[-3]], [nev, nev2], [1,3], [0, 205])
    ])
x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/1_percolation/1340_percolation.traj')
t0 = 68
plot(nbb[t0+1:]*1./x.get_nb()[t0:-1])
#plot(np.histogram(nev.values(), np.arange(t0+1, x.size+1))[0]*1./x.get_nb()[t0:-1])
plot(r_nev[t0:,0]-t0-1, r_nev[t0:,1])
#plot(np.histogram(
 #   [
  #      t for (a,b),t in nev.iteritems()
   #     if min(len(x.trajs[u]) for u in [a,b]) == x.size
    #    ], 
    #np.arange(t0+1, x.size+1))[0]*1./x.get_nb()[t0:-1])
plot(r_nevspan[t0:,0]-t0-1, r_nevspan[t0:,1])

for L in [5,13, 23]:
    rates = []
    for head,t0 in zip(['1_percolation/1340_percolation', '2_ageing/1415_ageing'], [68, 0]):
        x = xp.Experiment('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/%s.traj'%(head))
        lengths = [
            np.loadtxt(name, int) 
            for t, name in x.enum('_broken', 'length') 
            if t>0]
        maxlength = max([max(l) for l in lengths if len(l.shape)>0 and l.shape[0]>0])
        his = np.array([
            np.histogram(l, np.arange(-1, maxlength+1))[0]
            for l in lengths])
        rates.append(np.column_stack((
            np.arange(len(his)-t0),
            his[t0:,L:].sum(-1)*1./x.get_nb()[t0:-1]
            )))
    r = np.vstack((rates[0], rates[-1]*[3,1]+[205-68,0]))
    #10 points per decade
    li = np.unique(np.logspace(0, np.log10(len(r)), 10*np.log10(len(r))).astype(int))
    li = np.concatenate(([0],li))
    l = plot(r[li[:-1],0], [r[i:j,1].mean() for i,j in zip(li, li[1:])], label=lab)[0]



plot(np.arange(1,1000), 0.3*np.arange(1,1000)**(-1.5), 'k--')
xscale('log')
yscale('log')
plot([0]+range(2,10), [msd[10] for msd in msds])
plot([0]+range(2,10), [msd[10] for msd in msds])
xscale('log')
yscale('log')
plot([0]+range(2,10), [msd[10] for msd in msds])
xlabel(r'$\Lambda$')
ylabel(r'$\Delta r^2(\Delta t = 10\tau_B)/\sigma^2$')
#single trajectories that belong to a bond that never reform
nev_tr = np.unique(tr for bond in nev.iterkeys for tr in bond)
#index the trajectories belonging to a bond that never reform
tr2nev = -np.ones(len(x.trajs), int)
for itr, tr in enumerate(nev_tr):
    tr2nev[tr] = itr
#load positions
nev_p = np.zeros((x.size, len(nev_tr), 3))
pro = ProgressBar(x.size)
for t, name in x.enum():
    pos = np.loadtxt(name, skiprows=2)
    for itr, tr in enumerate(nev_tr):
        nev_p[t, itr] = pos[x.trajs[tr][t - x.starts[tr]]]
    pro.animate(t)
prebreak_bonds = np.array([len(nev), 3])
for bo, ((a,b), t) in enumerate(nev.iteritems):
    prebreak_bonds[bo] = np.diff(nev_p[t, tr2nev[[a,b]]], axis=0)
#single trajectories that belong to a bond that never reform
nev_tr = np.unique(tr for bond in nev.iterkeys() for tr in bond)
#index the trajectories belonging to a bond that never reform
tr2nev = -np.ones(len(x.trajs), int)
for itr, tr in enumerate(nev_tr):
    tr2nev[tr] = itr
#load positions
nev_p = np.zeros((x.size, len(nev_tr), 3))
pro = ProgressBar(x.size)
for t, name in x.enum():
    pos = np.loadtxt(name, skiprows=2)
    for itr, tr in enumerate(nev_tr):
        nev_p[t, itr] = pos[x.trajs[tr][t - x.starts[tr]]]
    pro.animate(t)
prebreak_bonds = np.array([len(nev), 3])
for bo, ((a,b), t) in enumerate(nev.iteritems):
    prebreak_bonds[bo] = np.diff(nev_p[t, tr2nev[[a,b]]], axis=0)
tr
tr, len(x.trajs)
x = xs[0]
#single trajectories that belong to a bond that never reform
nev_tr = np.unique(tr for bond in nev.iterkeys() for tr in bond)
#index the trajectories belonging to a bond that never reform
tr2nev = -np.ones(len(x.trajs), int)
for itr, tr in enumerate(nev_tr):
    tr2nev[tr] = itr
#load positions
nev_p = np.zeros((x.size, len(nev_tr), 3))
pro = ProgressBar(x.size)
for t, name in x.enum():
    pos = np.loadtxt(name, skiprows=2)
    for itr, tr in enumerate(nev_tr):
        nev_p[t, itr] = pos[x.trajs[tr][t - x.starts[tr]]]
    pro.animate(t)
prebreak_bonds = np.array([len(nev), 3])
for bo, ((a,b), t) in enumerate(nev.iteritems):
    prebreak_bonds[bo] = np.diff(nev_p[t, tr2nev[[a,b]]], axis=0)
nev_p.shape
pos.shape
len(x.trajs)
len(x.trajs[tr])
len(x.starts[tr])
x.starts[tr]
t
x = xs[0]
#single trajectories that belong to a bond that never reform
nev_tr = np.unique(tr for bond in nev.iterkeys() for tr in bond)
#index the trajectories belonging to a bond that never reform
tr2nev = -np.ones(len(x.trajs), int)
for itr, tr in enumerate(nev_tr):
    tr2nev[tr] = itr
#load positions
nev_p = np.zeros((x.size, len(nev_tr), 3))
pro = ProgressBar(x.size)
for t, name in x.enum():
    pos = np.loadtxt(name, skiprows=2)
    for itr, tr in enumerate(nev_tr):
        if t < x.starts[tr] or t - x.starts[tr] >= len(x.trajs[tr]):
            continue
        nev_p[t, itr] = pos[x.trajs[tr][t - x.starts[tr]]]
    pro.animate(t)
prebreak_bonds = np.array([len(nev), 3])
for bo, ((a,b), t) in enumerate(nev.iteritems):
    prebreak_bonds[bo] = np.diff(nev_p[t, tr2nev[[a,b]]], axis=0)
for bo, ((a,b), t) in enumerate(nev.iteritems()):
    prebreak_bonds[bo] = np.diff(nev_p[t, tr2nev[[a,b]]], axis=0)
prebreak_bonds[bo].shape
bo
prebreak_bonds = np.zeros([len(nev), 3])
for bo, ((a,b), t) in enumerate(nev.iteritems()):
    prebreak_bonds[bo] = np.diff(nev_p[t, tr2nev[[a,b]]], axis=0)
sta.plothist(np.sqrt(np.sum(prebreak_bonds**2)))
sta.plothist(np.sqrt(np.sum(prebreak_bonds**2, -1)))
prebreak_bonds = np.zeros([len(nev), 3])
for bo, ((a,b), t) in enumerate(nev.iteritems()):
    prebreak_bonds[bo] = np.diff(nev_p[t-1, tr2nev[[a,b]]], axis=0)
sta.plothist(np.sqrt(np.sum(prebreak_bonds**2, -1)))
#bonds just after breaking
postbreak_bonds = np.zeros([len(nev), 3])
for bo, ((a,b), t) in enumerate(nev.iteritems()):
    prebreak_bonds[bo] = np.diff(nev_p[t, tr2nev[[a,b]]], axis=0)
#bonds just before breaking
prebreak_bonds = np.zeros([len(nev), 3])
for bo, ((a,b), t) in enumerate(nev.iteritems()):
    prebreak_bonds[bo] = np.diff(nev_p[t-1, tr2nev[[a,b]]], axis=0)
#bonds just after breaking
postbreak_bonds = np.zeros([len(nev), 3])
for bo, ((a,b), t) in enumerate(nev.iteritems()):
    postbreak_bonds[bo] = np.diff(nev_p[t, tr2nev[[a,b]]], axis=0)
sta.plothist(np.sqrt(np.sum(prebreak_bonds**2, -1)))
sta.plothist(np.sqrt(np.sum(postbreak_bonds**2, -1)))
get_backend()
reload(sta)
sta.plothist(np.sqrt(np.sum(prebreak_bonds**2, -1)))
sta.plothist(np.sqrt(np.sum(postbreak_bonds**2, -1)))
sta.plothist(
get_ipython().magic(u'debug ')
sta.plothist(np.sqrt(np.sum(prebreak_bonds**2, -1)))
sta.plothist(np.sqrt(np.sum(postbreak_bonds**2, -1)))
sta.plothist(np.sqrt(np.sum(prebreak_bonds**2, -1)))
sta.plothist(np.sqrt(np.sum(postbreak_bonds**2, -1)))
sta.plothist(np.sqrt(np.dot(prebreak_bonds, postbreak_bonds)))
sta.plothist(np.sqrt(np.sum(prebreak_bonds**2, -1)))
sta.plothist(np.sqrt(np.sum(postbreak_bonds**2, -1)))
sta.plothist(np.sqrt(np.sum(prebreak_bonds * postbreak_bonds,-1)))
sta.plothist(np.sqrt(np.sum(prebreak_bonds**2, -1)))
sta.plothist(np.sqrt(np.sum(postbreak_bonds**2, -1)))
sta.plothist(np.sqrt(np.sum(prebreak_bonds / np.norm(prebreak_bonds) * postbreak_bonds,-1)))
sta.plothist(np.sqrt(np.sum(prebreak_bonds**2, -1)))
sta.plothist(np.sqrt(np.sum(postbreak_bonds**2, -1)))
sta.plothist(np.sqrt(np.sum(prebreak_bonds / np.sqrt(np.sum(prebreak_bonds**2)) * postbreak_bonds,-1)))
sta.plothist(np.sqrt(np.sum(prebreak_bonds**2, -1)))
sta.plothist(np.sqrt(np.sum(postbreak_bonds**2, -1)))
sta.plothist(np.sqrt(np.sum(prebreak_bonds * postbreak_bonds,-1))/ np.sqrt(np.sum(prebreak_bonds**2)*np.sum(postbreak_bonds**2)))
sta.plothist(np.sqrt(np.sum(prebreak_bonds * postbreak_bonds,-1))/ np.sqrt(np.sum(prebreak_bonds**2,-1)*np.sum(postbreak_bonds**2,-1)))
sta.plothist(np.sum(prebreak_bonds * postbreak_bonds,-1)/ np.sqrt(np.sum(prebreak_bonds**2,-1)*np.sum(postbreak_bonds**2,-1)))
for dt in range(10):
    normprod = []
    for bo, ((a,b), t) in enumerate(nev.iteritems()):
        if t + dt < min(len(x.trajs[u]) + x.starts[u] for u in (a,b)):
            v = np.diff(nev_p[t+dt, tr2nev[[a,b]]], axis=0)
            normprod.append(np.vdot(prebreak_bonds[bo], v) / np.sqrt(np.sum(prebreak_bonds[bo]**2) * np.vdot(v,v)))
    sta.plothist(np.arccos(normprod)/np.pi*180, np.arange(0,181))
for dt in range(10):
    normprod = []
    for bo, ((a,b), t) in enumerate(nev.iteritems()):
        if t + dt < min(len(x.trajs[u]) + x.starts[u] for u in (a,b)):
            v = np.diff(nev_p[t+dt, tr2nev[[a,b]]], axis=0)
            normprod.append(np.vdot(prebreak_bonds[bo], v) / np.sqrt(np.sum(prebreak_bonds[bo]**2) * np.vdot(v,v)))
    sta.plothist(np.arccos(normprod)/np.pi*180, np.arange(0,181), color=cm.autumn(dt/10.))
for dt in range(10):
    normprod = []
    for bo, ((a,b), t) in enumerate(nev.iteritems()):
        if t + dt < min(len(x.trajs[u]) + x.starts[u] for u in (a,b)):
            v = np.diff(nev_p[t+dt, tr2nev[[a,b]]], axis=0)
            normprod.append(np.vdot(prebreak_bonds[bo], v) / np.sqrt(np.sum(prebreak_bonds[bo]**2) * np.vdot(v,v)))
    h = np.histogram(np.arccos(normprod)/np.pi*180, np.arange(0,181))[0]
    plot(h, color=cm.autumn(dt/10.))
for dt in range(10):
    normprod = []
    for bo, ((a,b), t) in enumerate(nev.iteritems()):
        if t + dt < min(len(x.trajs[u]) + x.starts[u] for u in (a,b)):
            v = np.diff(nev_p[t+dt, tr2nev[[a,b]]], axis=0)
            normprod.append(np.vdot(prebreak_bonds[bo], v) / np.sqrt(np.sum(prebreak_bonds[bo]**2) * np.vdot(v,v)))
    h = np.histogram(np.arccos(normprod)/np.pi*180, np.arange(0,181))[0]
    plot(h/float(len(normprod)), color=cm.autumn(dt/10.))
for dt in 2**np.arange(5):
    normprod = []
    for bo, ((a,b), t) in enumerate(nev.iteritems()):
        if t + dt < min(len(x.trajs[u]) + x.starts[u] for u in (a,b)):
            v = np.diff(nev_p[t+dt, tr2nev[[a,b]]], axis=0)
            normprod.append(np.vdot(prebreak_bonds[bo], v) / np.sqrt(np.sum(prebreak_bonds[bo]**2) * np.vdot(v,v)))
    h = np.histogram(np.arccos(normprod)/np.pi*180, np.arange(0,181))[0]
    plot(h/float(len(normprod)), color=cm.autumn(dt/5.))
for dt in 2**np.arange(5):
    normprod = []
    for bo, ((a,b), t) in enumerate(nev.iteritems()):
        if t+dt < x.size/2 and t + dt < min(len(x.trajs[u]) + x.starts[u] for u in (a,b)):
            v = np.diff(nev_p[t+dt, tr2nev[[a,b]]], axis=0)
            normprod.append(np.vdot(prebreak_bonds[bo], v) / np.sqrt(np.sum(prebreak_bonds[bo]**2) * np.vdot(v,v)))
    h = np.histogram(np.arccos(normprod)/np.pi*180, np.arange(0,181))[0]
    plot(h/float(len(normprod)), color=cm.autumn(dt/5.))
for dt in 2**np.arange(5):
    normprod = []
    for bo, ((a,b), t) in enumerate(nev.iteritems()):
        if t < x.size/2 and t + dt < min(len(x.trajs[u]) + x.starts[u] for u in (a,b)):
            v = np.diff(nev_p[t+dt, tr2nev[[a,b]]], axis=0)
            normprod.append(np.vdot(prebreak_bonds[bo], v) / np.sqrt(np.sum(prebreak_bonds[bo]**2) * np.vdot(v,v)))
    h = np.histogram(np.arccos(normprod)/np.pi*180, np.arange(0,181))[0]
    plot(h/float(len(normprod)), color=cm.autumn(dt/5.))
for dt in 2**np.arange(5):
    normprod = []
    for bo, ((a,b), t) in enumerate(nev.iteritems()):
        if t < x.size/2 and t + dt < min(len(x.trajs[u]) + x.starts[u] for u in (a,b)):
            v = np.diff(nev_p[t+dt, tr2nev[[a,b]]], axis=0)
            normprod.append(np.vdot(prebreak_bonds[bo], v) / np.sqrt(np.sum(prebreak_bonds[bo]**2) * np.vdot(v,v)))
    h = np.histogram(np.arccos(normprod)/np.pi*180, np.arange(0,181))[0]
    plot(h/float(len(normprod)), color=cm.autumn(dt/5.))
xlabel('angle (degree)')
ylabel('probability')
co = np.zeros(x.size)
nb = np.zeros(x.size, int)
for bo, ((a,b), t) in enumerate(nev.iteritems()):
    tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
    vs = np.diff(nev_p[t:tmax, tr2nev[[a,b]]], axis=1)[:,0]
    co[:len(vs)] += np.sum(vs * prebreak_bonds[bo], -1)/ np.sqrt(np.sum(prebreak_bonds[bo]**2), np.sum(vs**2, -1))
    nb[:len(vs)] += 1
plot(co / np.maximum(1, nb))
tmax
t
vs.size
vs.shape
prebreak_bonds[bo].shape
(vs * prebreak_bonds[bo]).shape
np.sum(vs * prebreak_bonds[bo], -1).shape
co = np.zeros(x.size)
nb = np.zeros(x.size, int)
for bo, ((a,b), t) in enumerate(nev.iteritems()):
    tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
    vs = np.diff(nev_p[t:tmax, tr2nev[[a,b]]], axis=1)[:,0]
    co[:len(vs)] += np.sum(vs * prebreak_bonds[bo], -1)/ np.sqrt(np.sum(prebreak_bonds[bo]**2) * np.sum(vs**2, -1))
    nb[:len(vs)] += 1
plot(co / np.maximum(1, nb))
co = np.zeros(x.size)
nb = np.zeros(x.size, int)
for bo, ((a,b), t) in enumerate(nev.iteritems()):
    tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
    vs = np.diff(nev_p[t:tmax, tr2nev[[a,b]]], axis=1)[:,0]
    co[1:len(vs)+1] += np.sum(vs * prebreak_bonds[bo], -1)/ np.sqrt(np.sum(prebreak_bonds[bo]**2) * np.sum(vs**2, -1))
    nb[1:len(vs)+1] += 1
plot(co / np.maximum(1, nb))
co = np.zeros(x.size)
co[0] = 1
nb = np.zeros(x.size, int)
for bo, ((a,b), t) in enumerate(nev.iteritems()):
    tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
    vs = np.diff(nev_p[t:tmax, tr2nev[[a,b]]], axis=1)[:,0]
    co[1:len(vs)+1] += np.sum(vs * prebreak_bonds[bo], -1)/ np.sqrt(np.sum(prebreak_bonds[bo]**2) * np.sum(vs**2, -1))
    nb[1:len(vs)+1] += 1
plot(co / np.maximum(1, nb))
co = np.zeros(x.size)
co[0] = 1
nb = np.zeros(x.size, int)
for bo, ((a,b), t) in enumerate(nev.iteritems()):
    tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
    vs = np.diff(nev_p[t:tmax, tr2nev[[a,b]]], axis=1)[:,0]
    co[1:len(vs)+1] += np.sum(vs * prebreak_bonds[bo], -1)/ np.sqrt(np.sum(prebreak_bonds[bo]**2) * np.sum(vs**2, -1))
    nb[1:len(vs)+1] += 1
plot(co / np.maximum(1, nb))
xlabel(r'$\Delta t$')
xscale('log')
ylabel('vector correlation')
co = np.zeros(x.size)
co[0] = 1
nb = np.zeros(x.size, int)
for bo, ((a,b), t) in enumerate(nev.iteritems()):
    tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
    vs = np.diff(nev_p[t:tmax, tr2nev[[a,b]]], axis=1)[:,0]
    co[1:len(vs)+1] += np.sum(vs * prebreak_bonds[bo], -1)/ np.sqrt(np.sum(prebreak_bonds[bo]**2) * np.sum(vs**2, -1))
    nb[1:len(vs)+1] += 1
plot(co / np.maximum(1, nb))
xlabel(r'$\Delta t$')
xscale('log')
ylabel('vector correlation')
ylim(0,1)
co = np.zeros(x.size)
co[0] = 1
nb = np.zeros(x.size, int)
ts = [69,79,109,x.size]
for t1, t2 in zip(ts, ts[1:]):
    for bo, ((a,b), t) in enumerate(nev.iteritems()):
        if t<t1 or t>t2: continue
        tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
        vs = np.diff(nev_p[t:tmax, tr2nev[[a,b]]], axis=1)[:,0]
        co[1:len(vs)+1] += np.sum(vs * prebreak_bonds[bo], -1)/ np.sqrt(np.sum(prebreak_bonds[bo]**2) * np.sum(vs**2, -1))
        nb[1:len(vs)+1] += 1
    plot(co / np.maximum(1, nb))
xlabel(r'$\Delta t$')
xscale('log')
ylabel('vector correlation')
ylim(0,1)
co = np.zeros(x.size)
co[0] = 1
nb = np.zeros(x.size, int)
ts = [69,79,109,x.size]
for t1, t2 in zip(ts, ts[1:]):
    for bo, ((a,b), t) in enumerate(nev.iteritems()):
        if t<t1 or t>=t2: continue
        tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
        vs = np.diff(nev_p[t:tmax, tr2nev[[a,b]]], axis=1)[:,0]
        co[1:len(vs)+1] += np.sum(vs * prebreak_bonds[bo], -1)/ np.sqrt(np.sum(prebreak_bonds[bo]**2) * np.sum(vs**2, -1))
        nb[1:len(vs)+1] += 1
    plot(co / np.maximum(1, nb))
xlabel(r'$\Delta t$')
xscale('log')
ylabel('vector correlation')
ylim(0,1)
trajpos = [np.zeros([len(tr), 3]) for tr in x.trajs]
pro = ProgressBar(x.size)
for t, name in x.enum(ext='p2tr'):
    pos = np.loadtxt(x.get_format_string()%t, skiprows=2)
    p2tr = np.loadtxt(name, dtype=int)
    for p, tr in enumerate(p2tr):
        trajpos[tr][t - x.starts[tr]] = pos[tr]
    pro.animate(t)
trajpos = [np.zeros([len(tr), 3]) for tr in x.trajs]
pro = ProgressBar(x.size)
for t, name in x.enum(ext='p2tr'):
    pos = np.loadtxt(x.get_format_string()%t, skiprows=2)
    p2tr = np.loadtxt(name, dtype=int)
    for p, tr in enumerate(p2tr):
        trajpos[tr][t - x.starts[tr]] = pos[p]
    pro.animate(t)
nev_by_t = np.array([t, a, b] for (a,b), t in nev.iteritems())
nev_by_t = np.array([t, a, b] for (a,b), t in nev.iteritems())
nev_by_t = nev_by_t[np.argsort(nev_by_t[:,0])]
nev_by_t = np.array([[t, a, b] for (a,b), t in nev.iteritems()])
nev_by_t = nev_by_t[np.argsort(nev_by_t[:,0])]
nev_by_t = np.array([[t, a, b] for (a,b), t in nev.iteritems()])
nev_by_t = nev_by_t[np.argsort(nev_by_t[:,0])]
ti = [0] + list(np.where(np.diff(nev_by_t[:,0])>0)[0]) + [len(nev_by_t)]
nev_by_t = np.array([[t, a, b] for (a,b), t in nev.iteritems()])
nev_by_t = nev_by_t[np.argsort(nev_by_t[:,0])]
ti = [0] + list(1 + np.where(np.diff(nev_by_t[:,0])>0)[0]) + [len(nev_by_t)]
ti[:10]
ti[-10:]
nev_by_t = np.array([[t, a, b] for (a,b), t in nev.iteritems()])
nev_by_t = nev_by_t[np.argsort(nev_by_t[:,0])]
ix = [0] + list(1 + np.where(np.diff(nev_by_t[:,0])>0)[0]) + [len(nev_by_t)]
nev_by_t[ix[0]:ix[1], 0]
ngbbroke = [[] for t in range(x.size)]
for i, j in zip(ix, ix[1:]):
    t = nev_by_t[i,0] - 1 #time before breaking
    bonds = np.loadtxt(x.get_format_string(ext='bonds')%t, int)
    p2tr = np.loadtxt(x.get_format_string(ext='p2tr')%t, int)
    #construct neighbours for each particle
    ngb = [[] for n in range(len(p2tr))]
    for a,b in bonds:
        ngb[a].append(b)
        ngb[b].append(a)
    for tra, trb in nev_by_t[i:j,1:]: #iterate on broken bonds
        for tr in [tra, trb]:
            a = x.trajs[tr][t + x.starts[tr]]
            for ngb_tr in p2tr[ngb[a]]: #iterate on neighbours
                if ngb_tr not in [tra, trb]: #but not the broken one
                     ngbbroke[t+1].append(sorted([tr, ngb_tr])) #time after breaking, to be coherent with nev
ngbbroke = [[] for t in range(x.size)]
for i, j in zip(ix, ix[1:]):
    t = nev_by_t[i,0] - 1 #time before breaking
    bonds = np.loadtxt(x.get_format_string(ext='bonds')%t, int)
    p2tr = np.loadtxt(x.get_format_string(ext='p2tr')%t, int)
    #construct neighbours for each particle
    ngb = [[] for n in range(len(p2tr))]
    for a,b in bonds:
        ngb[a].append(b)
        ngb[b].append(a)
    for tra, trb in nev_by_t[i:j,1:]: #iterate on broken bonds
        for tr in [tra, trb]:
            a = x.trajs[tr][t - x.starts[tr]]
            for ngb_tr in p2tr[ngb[a]]: #iterate on neighbours
                if ngb_tr not in [tra, trb]: #but not the broken one
                     ngbbroke[t+1].append(sorted([tr, ngb_tr])) #time after breaking, to be coherent with nev
ngbbroke = [[] for t in range(x.size)]
pro = ProgressBar(x.size)
for i, j in zip(ix, ix[1:]):
    t = nev_by_t[i,0] - 1 #time before breaking
    bonds = np.loadtxt(x.get_format_string(ext='bonds')%t, int)
    p2tr = np.loadtxt(x.get_format_string(ext='p2tr')%t, int)
    #construct neighbours for each particle
    ngb = [[] for n in range(len(p2tr))]
    for a,b in bonds:
        ngb[a].append(b)
        ngb[b].append(a)
    for tra, trb in nev_by_t[i:j,1:]: #iterate on broken bonds
        for tr in [tra, trb]:
            a = x.trajs[tr][t - x.starts[tr]]
            for ngb_tr in p2tr[ngb[a]]: #iterate on neighbours
                if ngb_tr not in [tra, trb]: #but not the broken one
                     ngbbroke[t+1].append(sorted([tr, ngb_tr])) #time after breaking, to be coherent with nev
    pro.animate(t)
ts = [69,79,109,x.size]
for t1, t2 in zip(ts, ts[1:]):
    co = np.zeros(x.size)
    co[0] = 1
    nb = np.zeros(x.size, int)
    for bo, ((a,b), t) in enumerate(nev.iteritems()):
        if t<t1 or t>=t2: continue
        tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
        vs = np.diff(nev_p[t:tmax, tr2nev[[a,b]]], axis=1)[:,0]
        co[1:len(vs)+1] += np.sum(vs * prebreak_bonds[bo], -1)/ np.sqrt(np.sum(prebreak_bonds[bo]**2) * np.sum(vs**2, -1))
        nb[1:len(vs)+1] += 1
    plot(co / np.maximum(1, nb))
xlabel(r'$\Delta t$')
xscale('log')
ylabel('vector correlation')
ylim(0,1)
ts = [69,79,109,x.size]
for t1, t2 in zip(ts, ts[1:]):
    co = np.zeros(x.size)
    co[0] = 1
    nb = np.zeros(x.size, int)
    for bo, ((a,b), t) in enumerate(nev.iteritems()):
        if t<t1 or t>=t2: continue
        tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
        vs = np.diff(nev_p[t:tmax, tr2nev[[a,b]]], axis=1)[:,0]
        co[1:len(vs)+1] += np.sum(vs * prebreak_bonds[bo], -1)/ np.sqrt(np.sum(prebreak_bonds[bo]**2) * np.sum(vs**2, -1))
        nb[1:len(vs)+1] += 1
    plot(co / np.maximum(1, nb))
    

for t1, t2 in zip(ts, ts[1:]):
    co = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for t, bonds in enumerate(range(t1,t2), ngbbroke[t1:t2]):
        for a,b in bonds:
            tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
            vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
            vs /= np.sqrt(np.sum(vs**2, -1))[:, None]
            co[:len(vs)] += np.sum(vs * vs[0], -1)
            nb[:len(vs)] += 1
    plot(co / np.maximum(1, nb))
    
xlabel(r'$\Delta t$')
xscale('log')
ylabel('vector correlation')
ylim(0,1)
t1
t2
ts = [69,79,109,x.size]
for t1, t2 in zip(ts, ts[1:]):
    co = np.zeros(x.size)
    co[0] = 1
    nb = np.zeros(x.size, int)
    for bo, ((a,b), t) in enumerate(nev.iteritems()):
        if t<t1 or t>=t2: continue
        tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
        vs = np.diff(nev_p[t:tmax, tr2nev[[a,b]]], axis=1)[:,0]
        co[1:len(vs)+1] += np.sum(vs * prebreak_bonds[bo], -1)/ np.sqrt(np.sum(prebreak_bonds[bo]**2) * np.sum(vs**2, -1))
        nb[1:len(vs)+1] += 1
    plot(co / np.maximum(1, nb))
    

for t1, t2 in zip(ts, ts[1:]):
    co = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for t, bonds in zip(range(t1,t2), ngbbroke[t1:t2]):
        for a,b in bonds:
            tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
            vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
            vs /= np.sqrt(np.sum(vs**2, -1))[:, None]
            co[:len(vs)] += np.sum(vs * vs[0], -1)
            nb[:len(vs)] += 1
    plot(co / np.maximum(1, nb))
    
xlabel(r'$\Delta t$')
xscale('log')
ylabel('vector correlation')
ylim(0,1)
ts = [69,79,109,x.size]
for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    co = np.zeros(x.size)
    co[0] = 1
    nb = np.zeros(x.size, int)
    for bo, ((a,b), t) in enumerate(nev.iteritems()):
        if t<t1 or t>=t2: continue
        tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
        vs = np.diff(nev_p[t:tmax, tr2nev[[a,b]]], axis=1)[:,0]
        co[1:len(vs)+1] += np.sum(vs * prebreak_bonds[bo], -1)/ np.sqrt(np.sum(prebreak_bonds[bo]**2) * np.sum(vs**2, -1))
        nb[1:len(vs)+1] += 1
    plot(co / np.maximum(1, nb), c=c)
    

for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    co = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for t, bonds in zip(range(t1,t2), ngbbroke[t1:t2]):
        for a,b in bonds:
            tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
            vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
            vs /= np.sqrt(np.sum(vs**2, -1))[:, None]
            co[:len(vs)] += np.sum(vs * vs[0], -1)
            nb[:len(vs)] += 1
    plot(co / np.maximum(1, nb), ls=c+'--')
    
xlabel(r'$\Delta t$')
xscale('log')
ylabel('vector correlation')
ylim(0,1)
ts = [69,79,109,x.size]
for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    co = np.zeros(x.size)
    co[0] = 1
    nb = np.zeros(x.size, int)
    for bo, ((a,b), t) in enumerate(nev.iteritems()):
        if t<t1 or t>=t2: continue
        tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
        vs = np.diff(nev_p[t:tmax, tr2nev[[a,b]]], axis=1)[:,0]
        co[1:len(vs)+1] += np.sum(vs * prebreak_bonds[bo], -1)/ np.sqrt(np.sum(prebreak_bonds[bo]**2) * np.sum(vs**2, -1))
        nb[1:len(vs)+1] += 1
    plot(co / np.maximum(1, nb), c=c)
    

for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    co = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for t, bonds in zip(range(t1,t2), ngbbroke[t1:t2]):
        for a,b in bonds:
            tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
            vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
            vs /= np.sqrt(np.sum(vs**2, -1))[:, None]
            co[:len(vs)] += np.sum(vs * vs[0], -1)
            nb[:len(vs)] += 1
    plot(co / np.maximum(1, nb), ls='--', c=c)
    
xlabel(r'$\Delta t$')
xscale('log')
ylabel('vector correlation')
ylim(0,1)
ts = [69,79,109,x.size]
for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    co = np.zeros(x.size)
    co[0] = 1
    nb = np.zeros(x.size, int)
    for bo, ((a,b), t) in enumerate(nev.iteritems()):
        if t<t1 or t>=t2: continue
        tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
        vs = np.diff(nev_p[t:tmax, tr2nev[[a,b]]], axis=1)[:,0]
        co[1:len(vs)+1] += np.sum(vs * prebreak_bonds[bo], -1)/ np.sqrt(np.sum(prebreak_bonds[bo]**2) * np.sum(vs**2, -1))
        nb[1:len(vs)+1] += 1
    plot(co / np.maximum(1, nb), c=c)
    

for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    co = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for t, bonds in zip(range(t1,t2), ngbbroke[t1:t2]):
        for a,b in bonds:
            tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
            vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
            vs /= np.sqrt(np.sum(vs**2, -1))[:, None]
            co[:len(vs)] += np.sum(vs * vs[0], -1)
            nb[:len(vs)] += 1
    plot(co / np.maximum(1, nb), ls='--', c=c)
    
xlabel(r'$\Delta t$')
xscale('log')
ylabel('vector correlation')
ylim(0,1)
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/reorientation_postbreak_or_ngb.pdf')
ts = [69,79,109,x.size]
for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    sd = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for bo, ((a,b), t) in enumerate(nev.iteritems()):
        if t<t1 or t>=t2: continue
        tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
        vs = np.diff(nev_p[t-1:tmax, tr2nev[[a,b]]], axis=1)[:,0]
        sd[:len(vs)+1] += np.sum((vs - vs[0])**2, -1)
        sd[:len(vs)+1] += 1
    plot(sd / np.maximum(1, nb), c=c)
    

for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    sd = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for t, bonds in zip(range(t1,t2), ngbbroke[t1:t2]):
        for a,b in bonds:
            tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
            vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
            sd[:len(vs)] += np.sum((vs- vs[0])**2, -1)
            nb[:len(vs)] += 1
    plot(sd / np.maximum(1, nb), ls='--', c=c)
    
xlabel(r'$\Delta t$')
xscale('log')
yscale('log')
ylabel('msd')
ylim(0,1)
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/msd_postbreak_or_ngb.pdf')
ts = [69,79,109,x.size]
for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    sd = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for bo, ((a,b), t) in enumerate(nev.iteritems()):
        if t<t1 or t>=t2: continue
        tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
        vs = np.diff(nev_p[t-1:tmax, tr2nev[[a,b]]], axis=1)[:,0]
        sd[:len(vs)] += np.sum((vs - vs[0])**2, -1)
        sd[:len(vs)] += 1
    plot(sd / np.maximum(1, nb), c=c)
    

for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    sd = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for t, bonds in zip(range(t1,t2), ngbbroke[t1:t2]):
        for a,b in bonds:
            tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
            vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
            sd[:len(vs)] += np.sum((vs- vs[0])**2, -1)
            nb[:len(vs)] += 1
    plot(sd / np.maximum(1, nb), ls='--', c=c)
    
xlabel(r'$\Delta t$')
xscale('log')
yscale('log')
ylabel('msd')
ylim(0,1)
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/msd_postbreak_or_ngb.pdf')
ts = [69,79,109,x.size]
for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    sd = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for bo, ((a,b), t) in enumerate(nev.iteritems()):
        if t<t1 or t>=t2: continue
        tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
        vs = np.diff(nev_p[t-1:tmax, tr2nev[[a,b]]], axis=1)[:,0]
        sd[:len(vs)] += np.sum((vs - vs[0])**2, -1)
        sd[:len(vs)] += 1
    plot(sd / np.maximum(1, nb), c=c)
    

for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    sd = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for t, bonds in zip(range(t1,t2), ngbbroke[t1:t2]):
        for a,b in bonds:
            tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
            vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
            sd[:len(vs)] += np.sum((vs- vs[0])**2, -1)
            nb[:len(vs)] += 1
    plot(sd / np.maximum(1, nb), ls='--', c=c)
    
xlabel(r'$\Delta t$')
xscale('log')
yscale('log')
ylabel('msd')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/msd_postbreak_or_ngb.pdf')
ts = [69,79,109,x.size]
for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    sd = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for bo, ((a,b), t) in enumerate(nev.iteritems()):
        if t<t1 or t>=t2: continue
        tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
        vs = np.diff(nev_p[t-1:tmax, tr2nev[[a,b]]], axis=1)[:,0]
        sd[:len(vs)] += np.sum(vs**2, -1)
        sd[:len(vs)] += 1
    plot(sd / np.maximum(1, nb), c=c)
    

for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    sd = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for t, bonds in zip(range(t1,t2), ngbbroke[t1:t2]):
        for a,b in bonds:
            tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
            vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
            sd[:len(vs)] += np.sum(vs**2, -1)
            nb[:len(vs)] += 1
    plot(sd / np.maximum(1, nb), ls='--', c=c)
    
xlabel(r'$\Delta t$')
xscale('log')
yscale('log')
ylabel('msd')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/msd_postbreak_or_ngb.pdf')
ts = [69,79,109,x.size]
for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    sd = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for bo, ((a,b), t) in enumerate(nev.iteritems()):
        if t<t1 or t>=t2: continue
        tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
        vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
        sd[:len(vs)] += np.sum(vs**2, -1)
        nb[:len(vs)] += 1
    plot(sd / np.maximum(1, nb), c=c)
    

for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    sd = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for t, bonds in zip(range(t1,t2), ngbbroke[t1:t2]):
        for a,b in bonds:
            tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
            vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
            sd[:len(vs)] += np.sum(vs**2, -1)
            nb[:len(vs)] += 1
    plot(sd / np.maximum(1, nb), ls='--', c=c)
    
xlabel(r'$\Delta t$')
xscale('log')
yscale('log')
ylabel('msd')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/msd_postbreak_or_ngb.pdf')
ts = [69,79,109,x.size]
for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    sd = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for bo, ((a,b), t) in enumerate(nev.iteritems()):
        if t<t1 or t>=t2: continue
        tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
        vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
        sd[:len(vs)] += np.sum(vs**2, -1)
        nb[:len(vs)] += 1
    plot(sd / np.maximum(1, nb)/(2*x.radius)**2, c=c)
    

for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    sd = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for t, bonds in zip(range(t1,t2), ngbbroke[t1:t2]):
        for a,b in bonds:
            tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
            vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
            sd[:len(vs)] += np.sum(vs**2, -1)
            nb[:len(vs)] += 1
    plot(sd / np.maximum(1, nb)/(2*x.radius)**2, ls='--', c=c)
    
xlabel(r'$\Delta t$')
xscale('log')
yscale('log')
ylabel('msd')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/msd_postbreak_or_ngb.pdf')
ts = [69,79,109,x.size]
for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    sd = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for bo, ((a,b), t) in enumerate(nev.iteritems()):
        if t<t1 or t>=t2: continue
        tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
        vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
        sd[:len(vs)] += np.sum(vs**2, -1)
        nb[:len(vs)] += 1
    plot(sd / np.maximum(1, nb)/(2*x.rdf_radius())**2, c=c)
    

for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    sd = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for t, bonds in zip(range(t1,t2), ngbbroke[t1:t2]):
        for a,b in bonds:
            tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
            vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
            sd[:len(vs)] += np.sum(vs**2, -1)
            nb[:len(vs)] += 1
    plot(sd / np.maximum(1, nb)/(2*x.rdf_radius())**2, ls='--', c=c)
    
xlabel(r'$\Delta t$')
xscale('log')
yscale('log')
ylabel('msd')
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/msd_postbreak_or_ngb.pdf')
ts = [69,79,109,x.size]
for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    sd = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for bo, ((a,b), t) in enumerate(nev.iteritems()):
        if t<t1 or t>=t2: continue
        tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
        vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
        sd[:len(vs)] += np.sum(vs**2, -1)
        nb[:len(vs)] += 1
    plot(sd / np.maximum(1, nb)/(2*x.rdf_radius())**2, c=c)
    

for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    sd = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for t, bonds in zip(range(t1,t2), ngbbroke[t1:t2]):
        for a,b in bonds:
            tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
            vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
            sd[:len(vs)] += np.sum(vs**2, -1)
            nb[:len(vs)] += 1
    plot(sd / np.maximum(1, nb)/(2*x.rdf_radius())**2, ls='--', c=c)
    
xlabel(r'$\Delta t$')
xscale('log')
yscale('log')
ylabel('msd')
axhaxhline(1.3, 'k:')

savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/msd_postbreak_or_ngb.pdf')
ts = [69,79,109,x.size]
for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    sd = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for bo, ((a,b), t) in enumerate(nev.iteritems()):
        if t<t1 or t>=t2: continue
        tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
        vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
        sd[:len(vs)] += np.sum(vs**2, -1)
        nb[:len(vs)] += 1
    plot(sd / np.maximum(1, nb)/(2*x.rdf_radius())**2, c=c)
    

for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    sd = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for t, bonds in zip(range(t1,t2), ngbbroke[t1:t2]):
        for a,b in bonds:
            tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
            vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
            sd[:len(vs)] += np.sum(vs**2, -1)
            nb[:len(vs)] += 1
    plot(sd / np.maximum(1, nb)/(2*x.rdf_radius())**2, ls='--', c=c)
    
xlabel(r'$\Delta t$')
xscale('log')
yscale('log')
ylabel('msd')
axhline(1.3, 'k:')

savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/msd_postbreak_or_ngb.pdf')
get_ipython().magic(u'pinfo axhline')
ts = [69,79,109,x.size]
for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    sd = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for bo, ((a,b), t) in enumerate(nev.iteritems()):
        if t<t1 or t>=t2: continue
        tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
        vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
        sd[:len(vs)] += np.sum(vs**2, -1)
        nb[:len(vs)] += 1
    plot(sd / np.maximum(1, nb)/(2*x.rdf_radius())**2, c=c)
    

for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    sd = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for t, bonds in zip(range(t1,t2), ngbbroke[t1:t2]):
        for a,b in bonds:
            tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
            vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
            sd[:len(vs)] += np.sum(vs**2, -1)
            nb[:len(vs)] += 1
    plot(sd / np.maximum(1, nb)/(2*x.rdf_radius())**2, ls='--', c=c)
    
xlabel(r'$\Delta t$')
xscale('log')
yscale('log')
ylabel('msd')
axhline(1.3, c='k', ls=':')

savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/msd_postbreak_or_ngb.pdf')
ts = [69,79,109,x.size]
for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    sd = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for bo, ((a,b), t) in enumerate(nev.iteritems()):
        if t<t1 or t>=t2: continue
        tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
        vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
        sd[:len(vs)] += np.sum(vs**2, -1)
        nb[:len(vs)] += 1
    plot(sd / np.maximum(1, nb)/(2*x.rdf_radius())**2, c=c)
    

for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    sd = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for t, bonds in zip(range(t1,t2), ngbbroke[t1:t2]):
        for a,b in bonds:
            tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
            vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
            sd[:len(vs)] += np.sum(vs**2, -1)
            nb[:len(vs)] += 1
    plot(sd / np.maximum(1, nb)/(2*x.rdf_radius())**2, ls='--', c=c)
    
xlabel(r'$\Delta t$')
xscale('log')
yscale('log')
ylabel('msd')
axhline(1.3**2, c='k', ls=':')

savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/msd_postbreak_or_ngb.pdf')
ts = [69,79,109,x.size]
for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    sd = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for bo, ((a,b), t) in enumerate(nev.iteritems()):
        if t<t1 or t>=t2: continue
        tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
        vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
        sd[:len(vs)] += np.sum(vs**2, -1)
        nb[:len(vs)] += 1
    plot(sd / np.maximum(1, nb)/(2*x.rdf_radius())**2, c=c)
    

for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    sd = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for t, bonds in zip(range(t1,t2), ngbbroke[t1:t2]):
        for a,b in bonds:
            tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
            vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
            sd[:len(vs)] += np.sum(vs**2, -1)
            nb[:len(vs)] += 1
    plot(sd / np.maximum(1, nb)/(2*x.radius)**2, ls='--', c=c)
    
xlabel(r'$\Delta t$')
xscale('log')
yscale('log')
ylabel('msd')
axhline(1.3**2, c='k', ls=':')

savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/msd_postbreak_or_ngb.pdf')
ts = [69,79,109,x.size]
for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    sd = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for bo, ((a,b), t) in enumerate(nev.iteritems()):
        if t<t1 or t>=t2: continue
        tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
        vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
        sd[:len(vs)] += np.sum(vs**2, -1)
        nb[:len(vs)] += 1
    plot(sd / np.maximum(1, nb)/(2*x.rdf_radius())**2, c=c)
    

for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    sd = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for t, bonds in zip(range(t1,t2), ngbbroke[t1:t2]):
        for a,b in bonds:
            tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
            vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
            sd[:len(vs)] += np.sum(vs**2, -1)
            nb[:len(vs)] += 1
    plot(sd / np.maximum(1, nb)/(2*x.rdf_radius())**2, ls='--', c=c)
    
xlabel(r'$\Delta t$')
xscale('log')
yscale('log')
ylabel('msd')
axhline((1.3*x.radius/x.rdf_radius())**2, c='k', ls=':')

savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/msd_postbreak_or_ngb.pdf')
ts = [69,79,109,x.size]
for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    sd = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for bo, ((a,b), t) in enumerate(nev.iteritems()):
        if t<t1 or t>=t2: continue
        tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
        vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
        sd[:len(vs)] += np.sum(vs**2, -1)
        nb[:len(vs)] += 1
    plot(sd / np.maximum(1, nb)/(2*x.rdf_radius())**2, c=c)
    

for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    sd = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for t, bonds in zip(range(t1,t2), ngbbroke[t1:t2]):
        for a,b in bonds:
            tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
            vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
            sd[:len(vs)] += np.sum(vs**2, -1)
            nb[:len(vs)] += 1
    plot(sd / np.maximum(1, nb)/(2*x.rdf_radius())**2, ls='--', c=c)
    
xlabel(r'$\Delta t$')
xscale('log')
yscale('log')
ylabel('msd')
axhline((1.3*x.radius/x.rdf_radius())**2, c='k', ls=':')
tlog = np.logspace(0,3)
plot(tlog, 1.22**2*tlog**0.2, color='k', ls=':')

savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/msd_postbreak_or_ngb.pdf')
ts = [69,79,109,x.size]
for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    sd = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for bo, ((a,b), t) in enumerate(nev.iteritems()):
        if t<t1 or t>=t2: continue
        tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
        vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
        sd[:len(vs)] += np.sum(vs**2, -1)
        nb[:len(vs)] += 1
    plot(sd / np.maximum(1, nb)/(2*x.rdf_radius())**2, c=c)
    

for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    sd = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for t, bonds in zip(range(t1,t2), ngbbroke[t1:t2]):
        for a,b in bonds:
            tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
            vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
            sd[:len(vs)] += np.sum(vs**2, -1)
            nb[:len(vs)] += 1
    plot(sd / np.maximum(1, nb)/(2*x.rdf_radius())**2, ls='--', c=c)
    
xlabel(r'$\Delta t$')
xscale('log')
yscale('log')
ylabel('msd')
axhline((1.3*x.radius/x.rdf_radius())**2, c='k', ls=':')
tlog = np.logspace(0.5,1.5)
plot(tlog, 2**2*tlog**0.2, color='k', ls=':')

savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/msd_postbreak_or_ngb.pdf')
ts = [69,79,109,x.size]
for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    sd = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for bo, ((a,b), t) in enumerate(nev.iteritems()):
        if t<t1 or t>=t2: continue
        tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
        vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
        sd[:len(vs)] += np.sum(vs**2, -1)
        nb[:len(vs)] += 1
    plot(sd / np.maximum(1, nb)/(2*x.rdf_radius())**2, c=c)
    

for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    sd = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for t, bonds in zip(range(t1,t2), ngbbroke[t1:t2]):
        for a,b in bonds:
            tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
            vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
            sd[:len(vs)] += np.sum(vs**2, -1)
            nb[:len(vs)] += 1
    plot(sd / np.maximum(1, nb)/(2*x.rdf_radius())**2, ls='--', c=c)
    
xlabel(r'$\Delta t$')
xscale('log')
yscale('log')
ylabel('msd')
axhline((1.3*x.radius/x.rdf_radius())**2, c='k', ls=':')
tlog = np.logspace(0.5,1.5)
plot(tlog, 1.5**2*tlog**0.2, color='k', ls=':')

savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/msd_postbreak_or_ngb.pdf')
ts = [69,79,109,x.size]
for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    sd = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for bo, ((a,b), t) in enumerate(nev.iteritems()):
        if t<t1 or t>=t2: continue
        tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
        vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
        sd[:len(vs)] += np.sum(vs**2, -1)
        nb[:len(vs)] += 1
    plot(sd / np.maximum(1, nb)/(2*x.rdf_radius())**2, c=c)
    

for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    sd = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for t, bonds in zip(range(t1,t2), ngbbroke[t1:t2]):
        for a,b in bonds:
            tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
            vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
            sd[:len(vs)] += np.sum(vs**2, -1)
            nb[:len(vs)] += 1
    plot(sd / np.maximum(1, nb)/(2*x.rdf_radius())**2, ls='--', c=c)
    
xlabel(r'$\Delta t$')
xscale('log')
yscale('log')
ylabel('msd')
axhline((1.3*x.radius/x.rdf_radius())**2, c='k', ls=':')
tlog = np.logspace(0.5,1.5)
plot(tlog, 1.3**2*tlog**0.2, color='k', ls=':')

savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/msd_postbreak_or_ngb.pdf')
ts = [69,79,109,x.size]
for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    sd = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for bo, ((a,b), t) in enumerate(nev.iteritems()):
        if t<t1 or t>=t2: continue
        tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
        vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
        sd[:len(vs)] += np.sum(vs**2, -1)
        nb[:len(vs)] += 1
    plot(sd / np.maximum(1, nb)/(2*x.rdf_radius())**2, c=c)
    

for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    sd = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for t, bonds in zip(range(t1,t2), ngbbroke[t1:t2]):
        for a,b in bonds:
            tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
            vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
            sd[:len(vs)] += np.sum(vs**2, -1)
            nb[:len(vs)] += 1
    plot(sd / np.maximum(1, nb)/(2*x.rdf_radius())**2, ls='--', c=c)
    
xlabel(r'$\Delta t$')
xscale('log')
yscale('log')
ylabel('msd')
axhline((1.3*x.radius/x.rdf_radius())**2, c='k', ls=':')
tlog = np.logspace(0.5,1.5)
plot(tlog, 1.3**2*tlog**0.2, color='k', lw=2)

savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/msd_postbreak_or_ngb.pdf')
ts = [69,79,109,x.size]
for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    co = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for bo, ((a,b), t) in enumerate(nev.iteritems()):
        if t<t1 or t>=t2: continue
        tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
        #vs = np.diff(nev_p[t:tmax, tr2nev[[a,b]]], axis=1)[:,0]
        vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
        vs /= np.sqrt(np.sum(vs**2, -1))[:, None]
        #co[1:len(vs)+1] += np.sum(vs * prebreak_bonds[bo], -1)/ np.sqrt(np.sum(prebreak_bonds[bo]**2) * np.sum(vs**2, -1))
        co[:len(vs)] += np.sum(vs * vs[0], -1)
        nb[:len(vs)] += 1
    plot(co / np.maximum(1, nb), c=c)
    

for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    co = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for t, bonds in zip(range(t1,t2), ngbbroke[t1:t2]):
        for a,b in bonds:
            tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
            vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
            vs /= np.sqrt(np.sum(vs**2, -1))[:, None]
            co[:len(vs)] += np.sum(vs * vs[0], -1)
            nb[:len(vs)] += 1
    plot(co / np.maximum(1, nb), ls='--', c=c)
    
xlabel(r'$\Delta t$')
xscale('log')
ylabel('vector correlation')
ylim(0,1)
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/reorientation_postbreak_or_ngb.pdf')
ts = [69,79,109,x.size]
for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    co = np.zeros(x.size-t1)
    nb = np.zeros(x.size-t1, int)
    for bo, ((a,b), t) in enumerate(nev.iteritems()):
        if t<t1 or t>=t2: continue
        tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
        vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
        vs /= np.sqrt(np.sum(vs**2, -1))[:, None]
        co[:len(vs)] += np.sum(vs * vs[0], -1)
        nb[:len(vs)] += 1
    plot(co / np.maximum(1, nb), c=c)
    

for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    co = np.zeros(x.size-t1)
    nb = np.zeros(x.size-t1, int)
    for t, bonds in zip(range(t1,t2), ngbbroke[t1:t2]):
        for a,b in bonds:
            tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
            vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
            vs /= np.sqrt(np.sum(vs**2, -1))[:, None]
            co[:len(vs)] += np.sum(vs * vs[0], -1)
            nb[:len(vs)] += 1
    plot(co / np.maximum(1, nb), ls='--', c=c)
    
xlabel(r'$\Delta t$')
xscale('log')
ylabel('vector correlation')
ylim(0,1)
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/reorientation_postbreak_or_ngb.pdf')
ts = [69,79,109,x.size]
for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    co = np.zeros(x.size-t1+1)
    nb = np.zeros(x.size-t1+1, int)
    for bo, ((a,b), t) in enumerate(nev.iteritems()):
        if t<t1 or t>=t2: continue
        tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
        vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
        vs /= np.sqrt(np.sum(vs**2, -1))[:, None]
        co[:len(vs)] += np.sum(vs * vs[0], -1)
        nb[:len(vs)] += 1
    plot(co / np.maximum(1, nb), c=c)
    

for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    co = np.zeros(x.size)
    nb = np.zeros(x.size, int)
    for t, bonds in zip(range(t1,t2), ngbbroke[t1:t2]):
        for a,b in bonds:
            tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
            vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
            vs /= np.sqrt(np.sum(vs**2, -1))[:, None]
            co[:len(vs)] += np.sum(vs * vs[0], -1)
            nb[:len(vs)] += 1
    plot(co / np.maximum(1, nb), ls='--', c=c)
    
xlabel(r'$\Delta t$')
xscale('log')
ylabel('vector correlation')
ylim(0,1)
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/reorientation_postbreak_or_ngb.pdf')
ts = [69,79,109,x.size]
for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    co = np.zeros(x.size-t1+1)
    nb = np.zeros(x.size-t1+1, int)
    for bo, ((a,b), t) in enumerate(nev.iteritems()):
        if t<t1 or t>=t2: continue
        tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
        vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
        vs /= np.sqrt(np.sum(vs**2, -1))[:, None]
        co[:len(vs)] += np.sum(vs * vs[0], -1)
        nb[:len(vs)] += 1
    plot(co / np.maximum(1, nb), c=c)
    

for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    co = np.zeros(x.size-t1+1)
    nb = np.zeros(x.size-t1+1, int)
    for t, bonds in zip(range(t1,t2), ngbbroke[t1:t2]):
        for a,b in bonds:
            tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
            vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
            vs /= np.sqrt(np.sum(vs**2, -1))[:, None]
            co[:len(vs)] += np.sum(vs * vs[0], -1)
            nb[:len(vs)] += 1
    plot(co / np.maximum(1, nb), ls='--', c=c)
    
xlabel(r'$\Delta t$')
xscale('log')
ylabel('vector correlation')
ylim(0,1)
savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/reorientation_postbreak_or_ngb.pdf')
ts = [69,79,109,x.size]
for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    sd = np.zeros(x.size-t1+1)
    nb = np.zeros(x.size-t1+1, int)
    for bo, ((a,b), t) in enumerate(nev.iteritems()):
        if t<t1 or t>=t2: continue
        tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
        vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
        sd[:len(vs)] += np.sum(vs**2, -1)
        nb[:len(vs)] += 1
    plot(sd / np.maximum(1, nb)/(2*x.rdf_radius())**2, c=c)
    

for t1, t2, c in zip(ts, ts[1:], 'bgr'):
    sd = np.zeros(x.size-t1+1)
    nb = np.zeros(x.size-t1+1, int)
    for t, bonds in zip(range(t1,t2), ngbbroke[t1:t2]):
        for a,b in bonds:
            tmax = min(len(x.trajs[u]) + x.starts[u] for u in (a,b))
            vs = trajpos[a][t - 1 - x.starts[a]:tmax - x.starts[a]] - trajpos[b][t - 1 - x.starts[b]:tmax - x.starts[b]]
            sd[:len(vs)] += np.sum(vs**2, -1)
            nb[:len(vs)] += 1
    plot(sd / np.maximum(1, nb)/(2*x.rdf_radius())**2, ls='--', c=c)
    
xlabel(r'$\Delta t$')
xscale('log')
yscale('log')
ylabel('msd')
axhline((1.3*x.radius/x.rdf_radius())**2, c='k', ls=':')
tlog = np.logspace(0.5,1.5)
plot(tlog, 1.3**2*tlog**0.2, color='k', lw=2)

savefig('/media/storage/tsursawa/Prj1_Gelation/Thesis/0_Data/3_DenseGel/163A/msd_postbreak_or_ngb.pdf')
