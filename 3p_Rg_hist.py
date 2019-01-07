import numpy as np
from colloids import experiment as xp
from matplotlib import pyplot as plt

x = xp.Experiment('0_Data/0_Cluster/168A/168A_2004_percolation.traj')
Rgs = [
    np.sqrt(np.loadtxt(x.get_format_string()%t, skiprows=2)[
        np.array([list(map(int, line[:-1].split(','))) 
        for line in open(name) 
        if line.count(',')==2], int)
    ].var(axis=1).sum(axis=-1)) 
    for t, name in x.enum('_c2p', 'csv')
]

bins = np.linspace(0.5,1.02,51)
hists = np.array([np.histogram(Rg/x.radius/2, bins)[0] for Rg in Rgs])

fig = plt.figure('hist')
plt.clf()
#tree times close to t=0
for i in range(0,3):
    plt.plot(bins[:-1], hists[i*10:(i+1)*10].sum(0), c=plt.cm.autumn(i/2.), label='from %03d to %03d$\\tau_B$'%(i*10, (i+1)*10))
    np.savetxt('168A_3p_Rg_av10_t%03d.hist'%(10*i), np.column_stack((bins[:-1], hists[i*10:(i+1)*10].sum(0))), header='Rg number')

#long times
plt.plot(bins[:-1], hists[-150:].sum(0), c='k', label='from 350 to 500$\\tau_B$')
np.savetxt('168A_3p_Rg_av150_t%03d.hist'%(x.size-150), np.column_stack((bins[:-1], hists[-150:].sum(0))), header='Rg number')

plt.xlabel(r'$R_g/\sigma$')
plt.ylabel('#')
plt.savefig('3p_Rg_hist.pdf')
