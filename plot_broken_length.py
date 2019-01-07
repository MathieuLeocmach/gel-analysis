from matplotlib.pyplot import *
from scipy.optimize import leastsq
from colloids import statistics as sta

histlength150A = np.loadtxt('3_DenseGel/150A/2_ageing/broken_length.hist', usecols=[1])
histlength170A = np.loadtxt('2_MidGel/170A/2_ageing/broken_length.hist', usecols=[1])
histlength172A = np.loadtxt('0_Data/1_DiluteGel/172A/2_ageing/broken_length.hist', usecols=[1])

figure('broken length')
clf()

for hl, c, m,l in zip([histlength150A, histlength170A, histlength172A], 'bgr', 'vo^', ['dense', 'middle', 'dilute']):
    scatter(np.arange(2,len(hl)-1)-1, hl[3:], c=c, marker=m, label=l)
    h = sta.decimate_hist(np.column_stack((np.arange(2,len(hl)-1)-1, hl[3:])), col=1, thr=10)
    h = np.vstack((h, [len(hl)-2, h[-1,1:]]))
    step(h[:-1,0], h[:-1,1]/np.diff(h[:,0]).astype(float), where='post', c=c, label=l)
    
plot(2**np.arange(6), 1e6*(2**np.arange(6))**(-3.))
plot(2**np.arange(6), 1e4*(2**np.arange(6))**(-4.))

xscale('log');yscale('log')
xlim(0.9,140)
ylim(1,1e6)
ylabel('count')
xlabel('distance')
legend(loc='upper right')
savefig('broken_length_xylog.pdf')
