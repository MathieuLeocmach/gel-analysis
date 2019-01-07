from matplotlib.pyplot import *
from scipy.optimize import leastsq

fig=figure('Structure factor: cell vs capillary')
clf()

res = np.load('structure_factors/Res359A_scan2_structure_factor.npy')
cap = np.load('structure_factors/Cap359_scan2_structure_factor.npy')

for data, l, c, m in [[res, 'cell', 'b', 's'], [cap, 'capillary', 'r', 'o']]:
    q = np.arange(data.shape[1])*data[0,0]*10
    S = data.mean(0)
    scatter(q[3:], S[3:], label=l, c=c, marker=m)
    params = np.array([
        leastsq(
            lambda p, x, y : np.log(p[0]*x**(-p[1]))-np.log(y), 
            [0.6, 2.5], args=(q[qm:qM], S[qm:qM])
            )[0]
        for qm in range(9,12)
        for qM in range(29,31)
        ])
    plot(q[10:30], params[:,0].mean()*q[10:30]**(-params[:,1].mean()), color=c)
    print '%s: exponent = %.2f +- %.2f \t prefactor=%.2f +- %.2f'%(l, params[:,1].mean(), params[:,1].std(), params[:,0].mean(), params[:,0].std())

xscale('log');yscale('log')
xlabel(r'$q\sigma$')
ylabel(r'$S(q)$')
xlim(0.3, 20)
ylim(0.2,6)
legend(loc='lower left')
savefig('cell_vs_cap_Sq.pdf')
